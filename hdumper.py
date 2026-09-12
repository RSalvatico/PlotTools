import ROOT
import argparse
import glob
import csv
import os
from colorama import Fore, Style 
import numpy as np
import multiprocessing as mp
from functools import partial
import correctionlib
import flavTagWeightPlugin

ROOT.ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)
ROOT.TTreeCache.SetLearnEntries(200)
ROOT.gEnv.SetValue("TFile.AsyncPrefetching", 2)
ROOT.TH1.SetDefaultSumw2(True)

suffix_dict = {'base' : '', 'ttLF' : '_0', 'ttcj' : '_41', 'tt2c' : '_42', 'ttcc' : '_43', 'ttbj' : '_51', 'tt2b' : '_52', 'ttbb' : '_53'}
perProcessSysts = ["topHdampWeight_", "bFragWeight_", "bFragPetersonWeight_", "LHE_muF_", "LHE_muR_", "PS_fsr_", "PS_isr_", "minorBkg_PS_ISR_", "minorBkg_PS_FSR_"]


def add_overflow_underflow(hist):
    """
    Add underflow (bin 0) to first bin and overflow (bin N+1) to last bin.
    
    Parameters:
    - hist: ROOT.TH1 histogram
    """
    nbins = hist.GetNbinsX()
    
    # Add underflow to first bin
    underflow = hist.GetBinContent(0)
    underflow_err = hist.GetBinError(0)
    first_bin = hist.GetBinContent(1)
    first_bin_err = hist.GetBinError(1)
    
    hist.SetBinContent(1, first_bin + underflow)
    hist.SetBinError(1, np.sqrt(first_bin_err**2 + underflow_err**2))
    hist.SetBinContent(0, 0)
    hist.SetBinError(0, 0)
    
    # Add overflow to last bin
    overflow = hist.GetBinContent(nbins + 1)
    overflow_err = hist.GetBinError(nbins + 1)
    last_bin = hist.GetBinContent(nbins)
    last_bin_err = hist.GetBinError(nbins)
    
    hist.SetBinContent(nbins, last_bin + overflow)
    hist.SetBinError(nbins, np.sqrt(last_bin_err**2 + overflow_err**2))
    hist.SetBinContent(nbins + 1, 0)
    hist.SetBinError(nbins + 1, 0)
    
    return hist

def process_tree(infile, outfile, tree_name, hist_configs, year, selections, eventClassification, use5FS, count_events, do_systematics=False, flavtag_sf_json=None, flavtag_sf_name=None):
    """
    Processes a TTree, converts it to multiple TH1Ds for specified branches, and saves them to a ROOT file.

    Parameters:
    - input_files: List of input ROOT files.
    - output_files: List of output ROOT files.
    - tree_names: List of TTree names corresponding to input files.
    - hist_configs: List of dictionaries with keys 'branch', 'nbins', 'xmin', 'xmax'.
    - year: Year of data taking.
    - selections: String containing common event preselection.
    - eventClassification: Boolean indicating whether to apply event classification.
    - use5FS: Boolean indicating whether to use 5-flavor scheme MC for ttbb and ttbj processes.
    - count_events: Decide whether to count events for each selection (probably slows things down a bit).
    - flavtag_sf_json: Optional path to an alternate flavTaggingSF*.json.gz file; if given,
      flavTagWeight is recomputed on the fly from this file instead of the ntuple's branch.
    - flavtag_sf_name: Optional correctionlib correction name inside flavtag_sf_json.
    """
    print("")

    print(f"{Fore.RED}Processing file: {infile}{Style.RESET_ALL}")

    # Open input file
    input_file = ROOT.TFile.Open(infile)
    if not input_file or input_file.IsZombie():
        raise FileNotFoundError(f"Could not open file: {infile}")

    # Access the TTree
    tree = input_file.Get(tree_name)
    if not tree or not isinstance(tree, ROOT.TTree):
        raise ValueError(f"TTree '{tree_name}' not found in file '{infile}'.")
    
    # Optimize TTree reading
    tree.SetCacheSize(100000000)  # 100MB cache
    tree.AddBranchToCache("*", True)

    # Create RDataFrame from TTree
    df = ROOT.RDataFrame(tree)

    # Apply base selection everywhere (and early, to speed things up)
    base_filter = selections["base"]
    if "singlee" in infile:
        base_filter += " && passTrigMu==0" # Remove from the electron channel the events that fired the muon trigger. Could choose to do vice versa as well.
    df = df.Filter(base_filter)

    # Recompute the flavour-tagging weight (central value, plus every systematic variation
    # referenced by produce_systematics when running with --systematics) on the fly from an
    # alternate correctionlib SF file, instead of relying on the flavTagWeight* branches
    # already stored in the ntuple. The flavour-tagging entries do not depend on the
    # per-category suffix, so the needed variations are collected from the suffix-less dict.
    flavtag_weight_branch = flavTagWeightPlugin.STORED_WEIGHT_NAME
    if flavtag_sf_json and "data" not in infile.lower():
        needed_systs = flavTagWeightPlugin.extract_systematics(
            produce_systematics(year, '').values()) if do_systematics else None
        df, flavtag_weight_branch = flavTagWeightPlugin.define_flavtag_weights(
            df, year, json_path=flavtag_sf_json, correction_name=flavtag_sf_name,
            systematics=needed_systs)

    if eventClassification:
        print(f"{Fore.YELLOW}Running in event classification mode. Will define a series of fractional scores.{Style.RESET_ALL}")
        # Define the fractional scores
        df = df.Define("denominator", "score_ttbb + score_tt2b + score_ttbj + score_ttcc + score_tt2c + score_ttcj + score_ttLF") \
            .Define("fscore_ttbb", "score_ttbb / denominator") \
            .Define("fscore_tt2b", "score_tt2b / denominator") \
            .Define("fscore_ttbj", "score_ttbj / denominator") \
            .Define("fscore_ttcc", "score_ttcc / denominator") \
            .Define("fscore_tt2c", "score_tt2c / denominator") \
            .Define("fscore_ttcj", "score_ttcj / denominator") \
            .Define("fscore_ttLF", "score_ttLF / denominator")
    else:
        df = df.Define("ak4_1_pt", "ak4_pt.size() > 0 ? ak4_pt[0] : 0") \
            .Define("ak4_1_phi",   "ak4_phi.size() > 0 ? ak4_phi[0] : 0") \
            .Define("ak4_1_eta",   "ak4_eta.size() > 0 ? ak4_eta[0] : 0") \
            .Define("ak4_2_pt",    "ak4_pt.size() > 1 ? ak4_pt[1] : 0") \
            .Define("ak4_2_phi",   "ak4_phi.size() > 1 ? ak4_phi[1] : 0") \
            .Define("ak4_2_eta",   "ak4_eta.size() > 1 ? ak4_eta[1] : 0") \
            .Define("ak4_3_pt",    "ak4_pt.size() > 2 ? ak4_pt[2] : 0") \
            .Define("ak4_3_phi",   "ak4_phi.size() > 2 ? ak4_phi[2] : 0") \
            .Define("ak4_3_eta",   "ak4_eta.size() > 2 ? ak4_eta[2] : 0") \
            .Define("ak4_4_pt",    "ak4_pt.size() > 3 ? ak4_pt[3] : 0") \
            .Define("ak4_4_phi",   "ak4_phi.size() > 3 ? ak4_phi[3] : 0") \
            .Define("ak4_4_eta",   "ak4_eta.size() > 3 ? ak4_eta[3] : 0")

    tt_file_names = ["ttbb-4f", "ttbar-powheg"]
    tt4f_strings = ["ttbb", "ttbj", "tt2b"]
    tt_strings   = ["ttcc", "ttcj", "tt2c", "ttLF"]

    # Initialize counters for events
    local_total_MC_events = 0
    local_events_in_category = {key: 0 for key in selections.keys() if not key == "base"}

    flavTag_renormalization = {"base": "*1",
                               "ttbb": "*0.9393",
                               "tt2b": "*1.0024",
                               "ttbj": "*0.9541",
                               "ttcc": "*0.9908",
                               "tt2c": "*0.9388",
                               "ttcj": "*1.0094",
                               "ttLF": "*0.9974"
                               }

    histograms = {}
    # Process each selection-output combinations
    for selection_name in selections:

        # Apply base selection to every sample; apply the ttbar-specific selection to the right 4F, dps, and 5F powheg samples
        if not "base" in selection_name and not any(x in infile for x in tt_file_names): 
            continue
        if any(x in infile for x in tt_file_names) and "base" in selection_name:
            continue
        if any(x in selection_name for x in tt_strings) and not "powheg" in infile: 
            continue
        if use5FS: # only powheg
            if any(x in selection_name for x in tt4f_strings) and not "powheg" in infile:
                continue
        else: # dedicated 4fs samples for ttbb, ttbj, tt2b
            if any(x in selection_name for x in tt4f_strings) and not "bb" in infile:
                continue

        suffix = suffix_dict.get(selection_name, '')

        # Assign event weight based on data taking year and process type
        #weight = assign_event_weight(year, infile, suffix)
        systematics = produce_systematics(year, suffix)

        # Point the flavour-tagging terms at the recomputed columns, if we are using an
        # alternate SF file (no-op when flavtag_weight_branch is the stored branch name).
        systematics = {name: flavTagWeightPlugin.remap_expression(expr, flavtag_weight_branch)
                       for name, expr in systematics.items()}

        for syst in systematics.keys():
            if not do_systematics and not syst == "None":
                continue
            if syst == "None":
                weight = assign_event_weight(year, infile, suffix, flavtag_weight_branch=flavtag_weight_branch)
            else:
                weight = assign_event_weight(year, infile, suffix, systematics[syst], flavtag_weight_branch=flavtag_weight_branch)

            #if not "data" in infile and not "Data" in infile:
            #    weight = weight + flavTag_renormalization[selection_name] #FIXME
            #    print(f"weight is {weight}")
            if do_systematics:
                perProcessSystsWithoutLHEmuRmuF = [procDepSyst for procDepSyst in perProcessSysts if not (procDepSyst.startswith("LHE_muR") or procDepSyst.startswith("LHE_muF") or procDepSyst.startswith("minorBkg_PS_"))]
                if any(syst.startswith(procDepSyst) for procDepSyst in perProcessSystsWithoutLHEmuRmuF) and "tt" not in infile:
                    continue # Skip certain systematics if the process is not a ttbar one (including signal, ttH, and ttV)
                if syst.startswith("minorBkg") and "tt" in infile:
                    continue # Minor background systematics do not pertain to ttbar processes
    
            # Define a per-selection weight column to avoid re-defining the same column name
            weight_column = f"weight_{selection_name}_{syst}"
            if not "data" in infile and not "Data" in infile:
                print(f"Event weight: {weight}")
                #if "dps" in infile:
                #    weight = weight + "*4.52"
                df_weighted = df.Define(weight_column, weight)
            else:
                df_weighted = df.Define(weight_column, "1.")  # Set collision data weight to 1
    
            if not "base" in selection_name:
                print(f"Applying additional selection for {infile}: {Fore.RED}{selection_name}{Style.RESET_ALL}")
                ttbar_event_selection = f"{selections[selection_name]}"
                df_selected = df_weighted.Filter(ttbar_event_selection)
                if count_events:
                    print(f"Events passing additional ttbar selection: {df_selected.Count().GetValue()}")
            else:
                df_selected = df_weighted
    
            if not "Data" in infile and not "data" in infile and not "base" in selection_name:
                n_events = df_selected.Sum(weight_column).GetValue()
                local_total_MC_events += n_events
                local_events_in_category[selection_name] += n_events
    
            print(f"Applying selection: {Fore.GREEN}{selection_name}{Style.RESET_ALL}")
    
            # Define event classification for the dedicated mode
            if eventClassification:
                from configs.weights_and_constants import adhoc_selection, adhoc_binning
                adhoc_selection = adhoc_selection.copy()
                adhoc_binning = adhoc_binning.copy()
    
            # Create histograms for each branch
            final_df = dict()
            for hist_config in hist_configs:
                branch_name = hist_config['branch']
                nbins = int(hist_config['nbins'])
                xmin = float(hist_config['xmin'])
                xmax = float(hist_config['xmax'])
                print(f"Creating histogram for branch: {branch_name}")
                    
                final_df[branch_name] = df_selected.Filter(adhoc_selection[branch_name]) if eventClassification else df_selected
    
                hist_key = (branch_name, selection_name, syst)
                # Create histogram
                hist_name = f"h_{branch_name}_{syst}" if not syst == "None" else f"h_{branch_name}"
                if eventClassification:
                    #n_bins = 20 # Make many bins for these histograms. We will adjust them later.
                    histograms[hist_key] = final_df[branch_name].Histo1D((hist_name, f"Histogram of {branch_name}", len(adhoc_binning[branch_name])-1, adhoc_binning[branch_name]), branch_name, weight_column)
                    #histograms[hist_key] = final_df[branch_name].Histo1D((f"h_{branch_name}", f"Histogram of {branch_name}", nbins, xmin, xmax), branch_name, weight_column)
                else:
                    histograms[hist_key] = final_df[branch_name].Histo1D((hist_name, f"Histogram of {branch_name}", nbins, xmin, xmax), branch_name, weight_column)



    # Materialising histograms
    print(f"Materializing {len(histograms)} histograms...")
    materialized_hists = {}
    for key, hist_lazy in histograms.items():
        materialized_hists[key] = hist_lazy.GetPtr() 
        materialized_hists[key] = add_overflow_underflow(materialized_hists[key]) 

    output_file_handles = {}
    for key, hist in materialized_hists.items():
        branch_name, selection_name, syst = key
        tt_outfile_name = outfile.replace('.root','_'+selection_name+'.root')
        output_file = tt_outfile_name if not "base" in selection_name else outfile

        if output_file not in output_file_handles:
            output_file_handles[output_file] = ROOT.TFile(output_file, "RECREATE")

        output_file_handles[output_file].cd()
        hist.Write()
    
    # Close all files
    for fOut in output_file_handles.values():
        fOut.Close()

    input_file.Close()
    print(f"{Fore.GREEN}Completed processing {infile}{Style.RESET_ALL}")

    return (local_total_MC_events, local_events_in_category)


def process_trees_parallel(input_files, output_files, tree_name, hist_configs, year, selections, eventClassification, use5FS, count_events, do_systematics, flavtag_sf_json=None, flavtag_sf_name=None):
    """
    Basically a wrapper of process_tree to process multiple TTrees in parallel.
    """

    process_func = partial(
        process_tree,
        tree_name=tree_name,
        hist_configs=hist_configs,
        year=year,
        selections=selections,
        eventClassification=eventClassification,
        use5FS=use5FS,
        count_events=count_events,
        do_systematics=do_systematics,
        flavtag_sf_json=flavtag_sf_json,
        flavtag_sf_name=flavtag_sf_name
    )

    with mp.Pool(processes=min(len(input_files), mp.cpu_count())) as pool:
        results = pool.starmap(process_func, zip(input_files, output_files))

    # Aggregate results from all processes
    total_MC_events = 0
    events_in_category = {key: 0 for key in selections.keys() if not key == "base"}
    
    for local_total, local_category in results:
        total_MC_events += local_total
        for category, count in local_category.items():
            events_in_category[category] += count
    
    return (total_MC_events, events_in_category)

def read_csv(csv_file):
    """
    Open and read a csv file containing the name and the range of the variables to be plotted. 
    Fill in a list of dictionaries containing branch (i.e., variable name), nbins, xmin, and xmax information.

    Parameters:
    - csv_file: The csv file containing variable names and binning for the respective histograms.
    """
    with open(csv_file, mode = 'r') as f:
        csv_reader = csv.reader(f) 
        dict_list = [
            {'branch': line[0], 'nbins': line[1], 'xmin': line[2], 'xmax': line[3]}
            for line in csv_reader if not line[0] == 'Variable'
        ]
    # Note: the csv file must NOT contain empty lines.

    return dict_list

def assign_event_weight(year, infile, suffix, syst="", flavtag_weight_branch="flavTagWeight"):
    """
    Define the MC event weight according to the year. Collision data should be handled separately.

    Parameters:
    - year: Data taking year.
    - infile: Input file.
    - flavtag_weight_branch: Name of the branch/column to use for the flavour-tagging
      weight term (defaults to the "flavTagWeight" branch stored in the ntuple; pass
      the column added by flavTagWeightPlugin.define_flavtag_weight to use a weight
      recomputed on the fly from an alternate correctionlib SF file instead).
    """
    weight = "1"
    if year == 2024 or year == 2025:
        weight = f"lumiwgt*genWeight*xsecWeight*puWeight*muEffWeight*elEffWeight*{flavtag_weight_branch}*(((abs(lep1_pdgId)==11 && passTrigEl) || (abs(lep1_pdgId)==13 && passTrigMu)) && passmetfilters)"
    if "ttbar" in infile or "4f" in infile or "tt-vcb" in infile:
        weight = f"{weight}*TopPtWeight[1]*TopPtWeightNorm{suffix}[1]*TOPMLWeight[5]*TOPMLWeightNorm{suffix}[5]" # New variables in custom samples 

    if not syst == "":
        weight = f"{weight}*{syst}"
    
    return weight

def prepare_output(output_dir, input_files):
    """
    Prepare the output file names based on the input file names.

    Parameters:
    - output_dir: Output directory for the new ROOT files.
    - input_files: List of input ROOT files.
    """
    os.makedirs(output_dir, exist_ok=True)
    return [
        f"{output_dir}h_{input_file.split('/')[-1].replace('_tree.root','.root')}"
        for input_file in input_files
    ]

def merge_files(directory, input_files, output_file):
    """
    Merges multiple ROOT files into a single ROOT file.

    Parameters:
    - directory: Directory where the ROOT files are located.
    - input_files: List of input ROOT files.
    - output_file: Output ROOT file.
    """
    if not all([os.path.exists(directory+'/'+infile) for infile in input_files]):
        print(f"Input files {input_files} not found in directory: {directory}")
        return
    
    hadd_command = f"hadd -f {directory}/{output_file} {' '.join([directory+'/'+infile for infile in input_files])}"
    os.system(hadd_command)
    rm_command = f"rm {' '.join([directory+'/'+infile for infile in input_files])}"
    os.system(rm_command)

def produce_systematics(year, suffix):

    systematics = {"None" : "", 
                   #Pileup and lepton efficiencies
                   "CMS_pileup_%sUp"   % year  : "puWeightUp/puWeight", 
                   "CMS_pileup_%sDown" % year  : "puWeightDown/puWeight",
                   "CMS_trigEffUp"   : "trigEffWeightUp/trigEffWeight",
                   "CMS_trigEffDown" : "trigEffWeightDown/trigEffWeight",
                   "CMS_muEffUp"     : "muEffWeight_UP/muEffWeight",
                   "CMS_muEffDown"   : "muEffWeight_DOWN/muEffWeight",
                   "CMS_elEffUp"     : "elEffWeight_UP/elEffWeight",
                   "CMS_elEffDown"   : "elEffWeight_DOWN/elEffWeight",
                   "CMS_elSmearUp"   : "elSmear_UP",
                   "CMS_elSmearDown" : "elSmear_DOWN",
                   "CMS_elScaleUp"   : "elScale_UP",
                   "CMS_elScaleDown" : "elScale_DOWN",
                   "CMS_muSmearUp"   : "muSmear_UP",
                   "CMS_muSmearDown" : "muSmear_DOWN",
                   "CMS_muScaleUp"   : "muScale_UP",
                   "CMS_muScaleDown" : "muScale_DOWN",
                   # Flavor tagging
                   "CMS_flavTag_TTWeight_ttbarUp"     : "flavTagWeight_TTWeight_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_TTWeight_ttbarDown"   : "flavTagWeight_TTWeight_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_HDamp_ttbarUp"        : "flavTagWeight_HDamp_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_HDamp_ttbarDown"      : "flavTagWeight_HDamp_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_BDecay_ttbarUp"       : "flavTagWeight_BDecay_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_BDecay_ttbarDown"     : "flavTagWeight_BDecay_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_CDecay_ttbarUp"       : "flavTagWeight_CDecay_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_CDecay_ttbarDown"     : "flavTagWeight_CDecay_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_ttbarUp"         : "flavTagWeight_XSec_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_xsec_ttbarDown"       : "flavTagWeight_XSec_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_wjets_cUp"       : "flavTagWeight_XSec_WJets_c_UP/flavTagWeight",
                   "CMS_flavTag_xsec_wjets_cDown"     : "flavTagWeight_XSec_WJets_c_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_wjets_bUp"       : "flavTagWeight_XSec_WJets_b_UP/flavTagWeight",
                   "CMS_flavTag_xsec_wjets_bDown"     : "flavTagWeight_XSec_WJets_b_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_zjets_cUp"       : "flavTagWeight_XSec_ZJets_c_UP/flavTagWeight",
                   "CMS_flavTag_xsec_zjets_cDown"     : "flavTagWeight_XSec_ZJets_c_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_zjets_bUp"       : "flavTagWeight_XSec_ZJets_b_UP/flavTagWeight",
                   "CMS_flavTag_xsec_zjets_bDown"     : "flavTagWeight_XSec_ZJets_b_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_singlet_tChUp"   : "flavTagWeight_XSec_singlet_tCh_UP/flavTagWeight",
                   "CMS_flavTag_xsec_singlet_tChDown" : "flavTagWeight_XSec_singlet_tCh_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_singlet_tWUp"    : "flavTagWeight_XSec_singlet_tW_UP/flavTagWeight",
                   "CMS_flavTag_xsec_singlet_tWDown"  : "flavTagWeight_XSec_singlet_tW_DOWN/flavTagWeight",
                   "CMS_flavTag_xsec_VVUp"            : "flavTagWeight_XSec_VV_UP/flavTagWeight",
                   "CMS_flavTag_xsec_VVDown"          : "flavTagWeight_XSec_VV_DOWN/flavTagWeight",
                   "CMS_flavTag_PU_%sUp"     % year   : "flavTagWeight_PUWeight_UP/flavTagWeight",
                   "CMS_flavTag_PU_%sDown"   % year   : "flavTagWeight_PUWeight_DOWN/flavTagWeight",
                   "CMS_flavTag_Lumi_%sUp"   % year   : "flavTagWeight_Lumi_13p6TeV_%s_UP/flavTagWeight" % year,
                   "CMS_flavTag_Lumi_%sDown" % year   : "flavTagWeight_Lumi_13p6TeV_%s_DOWN/flavTagWeight" % year,
                   "CMS_flavTag_EleRecoUp"            : "flavTagWeight_Ele_Reco_UP/flavTagWeight",
                   "CMS_flavTag_EleRecoDown"          : "flavTagWeight_Ele_Reco_DOWN/flavTagWeight",
                   "CMS_flavTag_EleScaleUp"           : "flavTagWeight_Ele_Scale_UP/flavTagWeight",
                   "CMS_flavTag_EleScaleDown"         : "flavTagWeight_Ele_Scale_DOWN/flavTagWeight",
                   "CMS_flavTag_EleSmearUp"           : "flavTagWeight_Ele_Smear_UP/flavTagWeight",
                   "CMS_flavTag_EleSmearDown"         : "flavTagWeight_Ele_Smear_DOWN/flavTagWeight",
                   "CMS_flavTag_EleIDUp"              : "flavTagWeight_Ele_ID_UP/flavTagWeight",
                   "CMS_flavTag_EleIDDown"            : "flavTagWeight_Ele_ID_DOWN/flavTagWeight",
                   "CMS_flavTag_EleTriggerUp"         : "flavTagWeight_Ele_Trigger_UP/flavTagWeight",
                   "CMS_flavTag_EleTriggerDown"       : "flavTagWeight_Ele_Trigger_DOWN/flavTagWeight",
                   "CMS_flavTag_MuIDUp"               : "flavTagWeight_Mu_ID_UP/flavTagWeight",
                   "CMS_flavTag_MuIDDown"             : "flavTagWeight_Mu_ID_DOWN/flavTagWeight",
                   "CMS_flavTag_MuIsoUp"              : "flavTagWeight_Mu_Iso_UP/flavTagWeight",
                   "CMS_flavTag_MuIsoDown"            : "flavTagWeight_Mu_Iso_DOWN/flavTagWeight",
                   "CMS_flavTag_MuScaleUp"            : "flavTagWeight_Mu_Scale_UP/flavTagWeight",
                   "CMS_flavTag_MuScaleDown"          : "flavTagWeight_Mu_Scale_DOWN/flavTagWeight",
                   "CMS_flavTag_MuResolUp"            : "flavTagWeight_Mu_Resol_UP/flavTagWeight",
                   "CMS_flavTag_MuResolDown"          : "flavTagWeight_Mu_Resol_DOWN/flavTagWeight",
                   "CMS_flavTag_MuTriggerUp"          : "flavTagWeight_Mu_Trigger_UP/flavTagWeight",
                   "CMS_flavTag_MuTriggerDown"        : "flavTagWeight_Mu_Trigger_DOWN/flavTagWeight",
                   "CMS_flavTag_MET_UnclEnergyUp"     : "flavTagWeight_MET_UnclEnergy_UP/flavTagWeight",
                   "CMS_flavTag_MET_UnclEnergyDown"   : "flavTagWeight_MET_UnclEnergy_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C0_%sUp"   % year : "flavTagWeight_Stat_flavB_C0_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C0_%sDown" % year : "flavTagWeight_Stat_flavB_C0_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C1_%sUp"   % year : "flavTagWeight_Stat_flavB_C1_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C1_%sDown" % year : "flavTagWeight_Stat_flavB_C1_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C2_%sUp"   % year : "flavTagWeight_Stat_flavB_C2_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C2_%sDown" % year : "flavTagWeight_Stat_flavB_C2_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C3_%sUp"   % year : "flavTagWeight_Stat_flavB_C3_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C3_%sDown" % year : "flavTagWeight_Stat_flavB_C3_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C4_%sUp"   % year : "flavTagWeight_Stat_flavB_C4_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_C4_%sDown" % year : "flavTagWeight_Stat_flavB_C4_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B0_%sUp"   % year : "flavTagWeight_Stat_flavB_B0_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B0_%sDown" % year : "flavTagWeight_Stat_flavB_B0_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B1_%sUp"   % year : "flavTagWeight_Stat_flavB_B1_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B1_%sDown" % year : "flavTagWeight_Stat_flavB_B1_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B2_%sUp"   % year : "flavTagWeight_Stat_flavB_B2_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B2_%sDown" % year : "flavTagWeight_Stat_flavB_B2_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B3_%sUp"   % year : "flavTagWeight_Stat_flavB_B3_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B3_%sDown" % year : "flavTagWeight_Stat_flavB_B3_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B4_%sUp"   % year : "flavTagWeight_Stat_flavB_B4_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavB_B4_%sDown" % year : "flavTagWeight_Stat_flavB_B4_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C0_%sUp"   % year : "flavTagWeight_Stat_flavC_C0_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C0_%sDown" % year : "flavTagWeight_Stat_flavC_C0_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C1_%sUp"   % year : "flavTagWeight_Stat_flavC_C1_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C1_%sDown" % year : "flavTagWeight_Stat_flavC_C1_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C2_%sUp"   % year : "flavTagWeight_Stat_flavC_C2_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C2_%sDown" % year : "flavTagWeight_Stat_flavC_C2_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C3_%sUp"   % year : "flavTagWeight_Stat_flavC_C3_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C3_%sDown" % year : "flavTagWeight_Stat_flavC_C3_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C4_%sUp"   % year : "flavTagWeight_Stat_flavC_C4_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_C4_%sDown" % year : "flavTagWeight_Stat_flavC_C4_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B0_%sUp"   % year : "flavTagWeight_Stat_flavC_B0_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B0_%sDown" % year : "flavTagWeight_Stat_flavC_B0_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B1_%sUp"   % year : "flavTagWeight_Stat_flavC_B1_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B1_%sDown" % year : "flavTagWeight_Stat_flavC_B1_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B2_%sUp"   % year : "flavTagWeight_Stat_flavC_B2_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B2_%sDown" % year : "flavTagWeight_Stat_flavC_B2_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B3_%sUp"   % year : "flavTagWeight_Stat_flavC_B3_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B3_%sDown" % year : "flavTagWeight_Stat_flavC_B3_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B4_%sUp"   % year : "flavTagWeight_Stat_flavC_B4_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavC_B4_%sDown" % year : "flavTagWeight_Stat_flavC_B4_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C0_%sUp"   % year : "flavTagWeight_Stat_flavL_C0_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C0_%sDown" % year : "flavTagWeight_Stat_flavL_C0_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C1_%sUp"   % year : "flavTagWeight_Stat_flavL_C1_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C1_%sDown" % year : "flavTagWeight_Stat_flavL_C1_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C2_%sUp"   % year : "flavTagWeight_Stat_flavL_C2_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C2_%sDown" % year : "flavTagWeight_Stat_flavL_C2_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C3_%sUp"   % year : "flavTagWeight_Stat_flavL_C3_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C3_%sDown" % year : "flavTagWeight_Stat_flavL_C3_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C4_%sUp"   % year : "flavTagWeight_Stat_flavL_C4_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_C4_%sDown" % year : "flavTagWeight_Stat_flavL_C4_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B0_%sUp"   % year : "flavTagWeight_Stat_flavL_B0_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B0_%sDown" % year : "flavTagWeight_Stat_flavL_B0_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B1_%sUp"   % year : "flavTagWeight_Stat_flavL_B1_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B1_%sDown" % year : "flavTagWeight_Stat_flavL_B1_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B2_%sUp"   % year : "flavTagWeight_Stat_flavL_B2_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B2_%sDown" % year : "flavTagWeight_Stat_flavL_B2_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B3_%sUp"   % year : "flavTagWeight_Stat_flavL_B3_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B3_%sDown" % year : "flavTagWeight_Stat_flavL_B3_DOWN/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B4_%sUp"   % year : "flavTagWeight_Stat_flavL_B4_UP/flavTagWeight",
                   "CMS_flavTag_Stat_flavL_B4_%sDown" % year : "flavTagWeight_Stat_flavL_B4_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muF_ttbarUp"     : "flavTagWeight_LHEScaleWeight_muF_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muF_ttbarDown"   : "flavTagWeight_LHEScaleWeight_muF_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muR_ttbarUp"     : "flavTagWeight_LHEScaleWeight_muR_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muR_ttbarDown"   : "flavTagWeight_LHEScaleWeight_muR_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muF_singletUp"   : "flavTagWeight_LHEScaleWeight_muF_singlet_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muF_singletDown" : "flavTagWeight_LHEScaleWeight_muF_singlet_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muR_singletUp"   : "flavTagWeight_LHEScaleWeight_muR_singlet_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muR_singletDown" : "flavTagWeight_LHEScaleWeight_muR_singlet_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muF_wjetsUp"     : "flavTagWeight_LHEScaleWeight_muF_wjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muF_wjetsDown"   : "flavTagWeight_LHEScaleWeight_muF_wjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muR_wjetsUp"     : "flavTagWeight_LHEScaleWeight_muR_wjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muR_wjetsDown"   : "flavTagWeight_LHEScaleWeight_muR_wjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muF_zjetsUp"     : "flavTagWeight_LHEScaleWeight_muF_zjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muF_zjetsDown"   : "flavTagWeight_LHEScaleWeight_muF_zjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muR_zjetsUp"     : "flavTagWeight_LHEScaleWeight_muR_zjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muR_zjetsDown"   : "flavTagWeight_LHEScaleWeight_muR_zjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muF_dibosonUp"   : "flavTagWeight_LHEScaleWeight_muF_diboson_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muF_dibosonDown" : "flavTagWeight_LHEScaleWeight_muF_diboson_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_muR_dibosonUp"   : "flavTagWeight_LHEScaleWeight_muR_diboson_UP/flavTagWeight",
                   "CMS_flavTag_LHE_muR_dibosonDown" : "flavTagWeight_LHEScaleWeight_muR_diboson_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_ttbarUp"     : "flavTagWeight_LHEScaleWeight_PDF_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_ttbarDown"   : "flavTagWeight_LHEScaleWeight_PDF_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_singletUp"   : "flavTagWeight_LHEScaleWeight_PDF_singlet_UP/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_singletDown" : "flavTagWeight_LHEScaleWeight_PDF_singlet_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_wjetsUp"     : "flavTagWeight_LHEScaleWeight_PDF_wjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_wjetsDown"   : "flavTagWeight_LHEScaleWeight_PDF_wjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_zjetsUp"     : "flavTagWeight_LHEScaleWeight_PDF_zjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_zjetsDown"   : "flavTagWeight_LHEScaleWeight_PDF_zjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_dibosonUp"   : "flavTagWeight_LHEScaleWeight_PDF_diboson_UP/flavTagWeight",
                   "CMS_flavTag_LHE_PDF_dibosonDown" : "flavTagWeight_LHEScaleWeight_PDF_diboson_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_aS_ttbarUp"      : "flavTagWeight_LHEScaleWeight_aS_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_LHE_aS_ttbarDown"    : "flavTagWeight_LHEScaleWeight_aS_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_aS_singletUp"    : "flavTagWeight_LHEScaleWeight_aS_singlet_UP/flavTagWeight",
                   "CMS_flavTag_LHE_aS_singletDown"  : "flavTagWeight_LHEScaleWeight_aS_singlet_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_aS_wjetsUp"      : "flavTagWeight_LHEScaleWeight_aS_wjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_aS_wjetsDown"    : "flavTagWeight_LHEScaleWeight_aS_wjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_aS_zjetsUp"      : "flavTagWeight_LHEScaleWeight_aS_zjets_UP/flavTagWeight",
                   "CMS_flavTag_LHE_aS_zjetsDown"    : "flavTagWeight_LHEScaleWeight_aS_zjets_DOWN/flavTagWeight",
                   "CMS_flavTag_LHE_aS_dibosonUp"    : "flavTagWeight_LHEScaleWeight_aS_diboson_UP/flavTagWeight",
                   "CMS_flavTag_LHE_aS_dibosonDown"  : "flavTagWeight_LHEScaleWeight_aS_diboson_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_ISR_ttbarUp"      : "flavTagWeight_PSWeightISR_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_PS_ISR_ttbarDown"    : "flavTagWeight_PSWeightISR_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_FSR_ttbarUp"      : "flavTagWeight_PSWeightFSR_ttbar_UP/flavTagWeight",
                   "CMS_flavTag_PS_FSR_ttbarDown"    : "flavTagWeight_PSWeightFSR_ttbar_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_ISR_singletUp"    : "flavTagWeight_PSWeightISR_singlet_UP/flavTagWeight",
                   "CMS_flavTag_PS_ISR_singletDown"  : "flavTagWeight_PSWeightISR_singlet_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_FSR_singletUp"    : "flavTagWeight_PSWeightFSR_singlet_UP/flavTagWeight",
                   "CMS_flavTag_PS_FSR_singletDown"  : "flavTagWeight_PSWeightFSR_singlet_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_ISR_wjetsUp"      : "flavTagWeight_PSWeightISR_wjets_UP/flavTagWeight",
                   "CMS_flavTag_PS_ISR_wjetsDown"    : "flavTagWeight_PSWeightISR_wjets_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_FSR_wjetsUp"      : "flavTagWeight_PSWeightFSR_wjets_UP/flavTagWeight",
                   "CMS_flavTag_PS_FSR_wjetsDown"    : "flavTagWeight_PSWeightFSR_wjets_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_ISR_zjetsUp"      : "flavTagWeight_PSWeightISR_zjets_UP/flavTagWeight",
                   "CMS_flavTag_PS_ISR_zjetsDown"    : "flavTagWeight_PSWeightISR_zjets_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_FSR_zjetsUp"      : "flavTagWeight_PSWeightFSR_zjets_UP/flavTagWeight",
                   "CMS_flavTag_PS_FSR_zjetsDown"    : "flavTagWeight_PSWeightFSR_zjets_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_ISR_dibosonUp"    : "flavTagWeight_PSWeightISR_diboson_UP/flavTagWeight",
                   "CMS_flavTag_PS_ISR_dibosonDown"  : "flavTagWeight_PSWeightISR_diboson_DOWN/flavTagWeight",
                   "CMS_flavTag_PS_FSR_dibosonUp"    : "flavTagWeight_PSWeightFSR_diboson_UP/flavTagWeight",
                   "CMS_flavTag_PS_FSR_dibosonDown"  : "flavTagWeight_PSWeightFSR_diboson_DOWN/flavTagWeight",
                   "CMS_flavTag_JES_AbsoluteUp"      : "flavTagWeight_JESRegrouped_Absolute_UP/flavTagWeight",
                   "CMS_flavTag_JES_AbsoluteDown"    : "flavTagWeight_JESRegrouped_Absolute_DOWN/flavTagWeight",
                   "CMS_flavTag_JES_BBEC1Up"         : "flavTagWeight_JESRegrouped_BBEC1_UP/flavTagWeight",
                   "CMS_flavTag_JES_BBEC1Down"       : "flavTagWeight_JESRegrouped_BBEC1_DOWN/flavTagWeight",
                   "CMS_flavTag_JES_FlavorQCDUp"     : "flavTagWeight_JESRegrouped_FlavorQCD_UP/flavTagWeight",
                   "CMS_flavTag_JES_FlavorQCDDown"   : "flavTagWeight_JESRegrouped_FlavorQCD_DOWN/flavTagWeight",
                   "CMS_flavTag_JES_RelativeBalUp"   : "flavTagWeight_JESRegrouped_RelativeBal_UP/flavTagWeight",
                   "CMS_flavTag_JES_RelativeBalDown" : "flavTagWeight_JESRegrouped_RelativeBal_DOWN/flavTagWeight",
                   "CMS_flavTag_JES_Absolute_%sUp"   % year : "flavTagWeight_JESRegrouped_Absolute_%s_UP/flavTagWeight" % year,
                   "CMS_flavTag_JES_Absolute_%sDown" % year : "flavTagWeight_JESRegrouped_Absolute_%s_DOWN/flavTagWeight" % year,
                   "CMS_flavTag_JES_BBEC1_%sUp"      % year : "flavTagWeight_JESRegrouped_BBEC1_%s_UP/flavTagWeight" % year,
                   "CMS_flavTag_JES_BBEC1_%sDown"    % year : "flavTagWeight_JESRegrouped_BBEC1_%s_DOWN/flavTagWeight" % year,
                   "CMS_flavTag_JES_RelativeSample_%sUp"   % year : "flavTagWeight_JESRegrouped_RelativeSample_%s_UP/flavTagWeight" % year,
                   "CMS_flavTag_JES_RelativeSample_%sDown" % year : "flavTagWeight_JESRegrouped_RelativeSample_%s_DOWN/flavTagWeight" % year,
                   "CMS_flavTag_JEReta0to1p93_%sUp"   % year : "flavTagWeight_JEReta0to1p93_UP/flavTagWeight",
                   "CMS_flavTag_JEReta0to1p93_%sDown" % year : "flavTagWeight_JEReta0to1p93_DOWN/flavTagWeight",
                   "CMS_flavTag_JEReta1p93to2p5_%sUp"   % year : "flavTagWeight_JEReta1p93to2p5_UP/flavTagWeight",
                   "CMS_flavTag_JEReta1p93to2p5_%sDown" % year : "flavTagWeight_JEReta1p93to2p5_DOWN/flavTagWeight",
                   # Hdamp, b fragmentation, LHE scale, PS weights
                   "topHdampWeight_%sUp"   % year : f"TOPMLWeight[1]*TOPMLWeightNorm{suffix}[1]",
                   "topHdampWeight_%sDown" % year : f"TOPMLWeight[3]*TOPMLWeightNorm{suffix}[3]",
                   "bFragWeight_%sUp"   % year : f"(TOPMLWeight[4]*TOPMLWeightNorm{suffix}[4])/(TOPMLWeight[5]*TOPMLWeightNorm{suffix}[5])", #Divide by bFrag nominal and multiply by bFrag up
                   "bFragWeight_%sDown" % year : f"1/(TOPMLWeight[5]*TOPMLWeightNorm{suffix}[5])", #The standard samples are effectively bFrag down, so here just dividing by bFrag nominal and its renorm weight
                   "bFragPetersonWeight_%sUp"   % year : f"(bFragAndDecayWeight[3]*BFragAndDecayWeightNorm{suffix}[3])/(TOPMLWeight[5]*TOPMLWeightNorm{suffix}[5])", #A one-sided systematic
                   "bFragPetersonWeight_%sDown" % year : f"1.", #Effectively a one-sided systematic
                   # LHE for minor bkgs
                   #"LHE_minorBkg_muF_%sUp" % year : "LHEScaleWeight[5]*LHEScaleWeightNorm[5]",
                   #"LHE_minorBkg_muF_%sDown" % year : "LHEScaleWeight[3]*LHEScaleWeightNorm[3]",
                   #"LHE_minorBkg_muR_%sUp" % year : "LHEScaleWeight[7]*LHEScaleWeightNorm[7]",
                   #"LHE_minorBkg_muR_%sDown" % year : "LHEScaleWeight[1]*LHEScaleWeightNorm[1]",
                   # LHE for large bkgs
                   "LHE_muF_%sUp"   % year : f"LHEScaleWeight[5]*LHEScaleWeightNorm{suffix}[5]",
                   "LHE_muF_%sDown" % year : f"LHEScaleWeight[3]*LHEScaleWeightNorm{suffix}[3]",
                   "LHE_muR_%sUp"   % year : f"LHEScaleWeight[7]*LHEScaleWeightNorm{suffix}[7]",
                   "LHE_muR_%sDown" % year : f"LHEScaleWeight[1]*LHEScaleWeightNorm{suffix}[1]",
                   # PS for minor bkgs
                   "minorBkg_PS_ISR_%sUp"   % year : f"PSWeight[0]*PSWeightNorm{suffix}[0]",
                   "minorBkg_PS_ISR_%sDown" % year : f"PSWeight[2]*PSWeightNorm{suffix}[2]",
                   "minorBkg_PS_FSR_%sUp"   % year : f"PSWeight[1]*PSWeightNorm{suffix}[1]",
                   "minorBkg_PS_FSR_%sDown" % year : f"PSWeight[3]*PSWeightNorm{suffix}[3]",
                   # PS-fsr for large bkgs
                   "PS_fsr_G2GG_muR_%sDown" %year : f"PSWeight[6]*PSWeightNorm{suffix}[6]",
                   "PS_fsr_G2GG_muR_%sUp"   %year : f"PSWeight[7]*PSWeightNorm{suffix}[7]",
                   "PS_fsr_G2QQ_muR_%sDown" %year : f"PSWeight[8]*PSWeightNorm{suffix}[8]",
                   "PS_fsr_G2QQ_muR_%sUp"   %year : f"PSWeight[9]*PSWeightNorm{suffix}[9]",                              
                   "PS_fsr_Q2QG_muR_%sDown" %year : f"PSWeight[10]*PSWeightNorm{suffix}[10]",
                   "PS_fsr_Q2QG_muR_%sUp"   %year : f"PSWeight[11]*PSWeightNorm{suffix}[11]",                              
                   "PS_fsr_X2XG_muR_%sDown" %year : f"PSWeight[12]*PSWeightNorm{suffix}[12]",
                   "PS_fsr_X2XG_muR_%sUp"   %year : f"PSWeight[13]*PSWeightNorm{suffix}[13]",
                   "PS_fsr_G2GG_cNS_%sDown" %year : f"PSWeight[14]*PSWeightNorm{suffix}[14]",
                   "PS_fsr_G2GG_cNS_%sUp"   %year : f"PSWeight[15]*PSWeightNorm{suffix}[15]",
                   "PS_fsr_G2QQ_cNS_%sDown" %year : f"PSWeight[16]*PSWeightNorm{suffix}[16]",
                   "PS_fsr_G2QQ_cNS_%sUp"   %year : f"PSWeight[17]*PSWeightNorm{suffix}[17]",
                   "PS_fsr_G2QG_cNS_%sDown" %year : f"PSWeight[18]*PSWeightNorm{suffix}[18]",
                   "PS_fsr_G2QG_cNS_%sUp"   %year : f"PSWeight[19]*PSWeightNorm{suffix}[19]",
                   "PS_fsr_X2XG_cNS_%sDown" %year : f"PSWeight[20]*PSWeightNorm{suffix}[20]",
                   "PS_fsr_X2XG_cNS_%sUp"   %year : f"PSWeight[21]*PSWeightNorm{suffix}[21]",
                   # PS-isr for large bkgs
                   "PS_isr_G2GG_muR_%sDown" %year : f"PSWeight[28]*PSWeightNorm{suffix}[28]",
                   "PS_isr_G2GG_muR_%sUp"   %year : f"PSWeight[29]*PSWeightNorm{suffix}[29]",
                   "PS_isr_G2QQ_muR_%sDown" %year : f"PSWeight[30]*PSWeightNorm{suffix}[30]",
                   "PS_isr_G2QQ_muR_%sUp"   %year : f"PSWeight[31]*PSWeightNorm{suffix}[31]",
                   "PS_isr_Q2QG_muR_%sDown" %year : f"PSWeight[32]*PSWeightNorm{suffix}[32]",
                   "PS_isr_Q2QG_muR_%sUp"   %year : f"PSWeight[33]*PSWeightNorm{suffix}[33]",
                   "PS_isr_X2XG_muR_%sDown" %year : f"PSWeight[34]*PSWeightNorm{suffix}[34]",
                   "PS_isr_X2XG_muR_%sUp"   %year : f"PSWeight[35]*PSWeightNorm{suffix}[35]",
                   "PS_isr_G2GG_cNS_%sDown" %year : f"PSWeight[36]*PSWeightNorm{suffix}[36]",
                   "PS_isr_G2GG_cNS_%sUp"   %year : f"PSWeight[37]*PSWeightNorm{suffix}[37]",
                   "PS_isr_G2QQ_cNS_%sDown" %year : f"PSWeight[38]*PSWeightNorm{suffix}[38]",
                   "PS_isr_G2QQ_cNS_%sUp"   %year : f"PSWeight[39]*PSWeightNorm{suffix}[39]",
                   "PS_isr_G2QG_cNS_%sDown" %year : f"PSWeight[40]*PSWeightNorm{suffix}[40]",
                   "PS_isr_G2QG_cNS_%sUp"   %year : f"PSWeight[41]*PSWeightNorm{suffix}[41]",
                   "PS_isr_X2XG_cNS_%sDown" %year : f"PSWeight[42]*PSWeightNorm{suffix}[42]",
                   "PS_isr_X2XG_cNS_%sUp"   %year : f"PSWeight[43]*PSWeightNorm{suffix}[43]",
               }
               
    return systematics
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process ROOT TTrees into TH1D histograms.")
    parser.add_argument("--input_dirs", nargs='+', required=True, help="List of directories where the ROOT files are fetched.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory of the new ROOT files.")
    parser.add_argument("--tree_name", type=str, required=True, help="List of TTree names in the input files.")
    parser.add_argument("--nbins", type=int, required=False, help="Number of bins for the histograms.")
    parser.add_argument("--xmin", type=float, required=False, help="Minimum value for the histograms.")
    parser.add_argument("--xmax", type=float, required=False, help="Maximum value for the histograms.")
    parser.add_argument("--input_csv", type=str, required=True, help="The csv file to read variables and ranges from.")
    parser.add_argument("--year", type=int, required=True, help="Data taking year.")
    parser.add_argument("--electron", nargs="?", const=1, type=bool, default=False, required=False, help="Process electron channel only.")
    parser.add_argument("--muon", nargs="?", const=1, type=bool, default=False, required=False, help="Process muon channel only.")
    parser.add_argument("--add_selection", type=str, required=False, help="Additional selection to apply to all processes.")
    parser.add_argument("--count_events", nargs="?", const=1, type=bool, default=False, required=False, help="Count events for each selection.")
    parser.add_argument("--eventClassification", nargs="?", const=1, type=bool, default=False, required=False, help="Apply event classification selection.")
    parser.add_argument("--use5FS", nargs="?", const=1, type=bool, default=False, required=False, help="Use 5-flavor scheme.")
    parser.add_argument("--systematics", nargs="?", const=1, type=bool, default=False, required=False, help="Make systematic variations.")
    parser.add_argument("--flavtag_sf_json", type=str, required=False, default=None,
                         help="Path to an alternate flavTaggingSF*.json.gz correctionlib file. "
                              "If given, flavTagWeight is recomputed on the fly from this file "
                              "instead of using the flavTagWeight branch stored in the ntuple.")
    parser.add_argument("--flavtag_sf_name", type=str, required=False, default=None,
                         help="Name of the correction inside --flavtag_sf_json "
                              "(defaults to the standard per-year name, e.g. particleNetAK4_shape "
                              "or UParTAK4_pseudocontinuous).")


    args = parser.parse_args()

    # Get input files from the input_dirs list
    input_files = []
    for input_dir in args.input_dirs:
        input_files.extend(glob.glob(f"{input_dir}*.root"))
    input_files = sorted(set(input_files))

    # Prepare list of output files based on the name of the input files
    output_files = prepare_output(args.output_dir, input_files)

    # Prepare histogram configurations for each branch
    hist_configs = read_csv(args.input_csv)

    selections = {"base": "n_ak4>=4 && n_btagM>=2 && n_ctagM>=1",
                 "ttbb" : "genEventClassifier==9",
                 "ttbj" : "genEventClassifier==7",
                 "tt2b" : "genEventClassifier==8",
                 "ttcc" : "genEventClassifier==6",
                 "ttcj" : "genEventClassifier==4",
                 "tt2c" : "genEventClassifier==5",
                 "ttLF" : "tt_category==0"
    }

    # Apply trigger selection to separate channels if requested
    if args.electron:
        selections["base"] += " && passTrigEl"
    if args.muon:
        selections["base"] += " && passTrigMu"

    use5FS = False
    if args.use5FS:
        use5FS = True
        print(f"{Fore.GREEN}Using 5-flavor scheme for ttbb, tt2b, and ttbj processes.{Style.RESET_ALL}")

    # Apply additional selections if specified
    if args.add_selection:
        for key in selections.keys():
            selections[key] += f" && ({args.add_selection})"

    print(f"{Fore.YELLOW}Final selections to be applied:{Style.RESET_ALL}")
    for key, value in selections.items():
        print(f"{Fore.YELLOW} - {key}: {value}{Style.RESET_ALL}")

    # Process the trees and get event counts
    total_MC_events, events_in_category = process_trees_parallel(input_files, output_files, args.tree_name, hist_configs, args.year, selections, args.eventClassification, use5FS, args.count_events, args.systematics, args.flavtag_sf_json, args.flavtag_sf_name)

    # Make sure the previous step is completed before merging files
    ROOT.gSystem.Exec("sync")

    # Merge some of the output files
    ttH_list = ["h_ttHbb.root", "h_ttHcc.root", "h_ttZ.root", "h_ttW.root", "h_diboson.root", "h_singletop.root", "h_wjets.root"]
    merge_files(args.output_dir, ttH_list, "h_others.root")

    merge_files(args.output_dir, ["h_ttbb-4f_ttbb.root", "h_ttbb-4f-dps_ttbb.root"], "h_ttbb.root")
    merge_files(args.output_dir, ["h_ttbb-4f_tt2b.root", "h_ttbb-4f-dps_tt2b.root"], "h_tt2b.root")
    merge_files(args.output_dir, ["h_ttbb-4f_ttbj.root", "h_ttbb-4f-dps_ttbj.root"], "h_ttbj.root")

    data_list = ["h_singlee.root", "h_singlemu.root"]
    merge_files(args.output_dir, data_list, "h_Data.root")

    if args.count_events:
        print(f"{Fore.YELLOW}Total MC events after preselection: {total_MC_events}{Style.RESET_ALL}")
        print(f"{Fore.YELLOW}Event counts in each category:{Style.RESET_ALL}")
        for category, count in events_in_category.items():
            print(f"{Fore.YELLOW} - {category}: {count}{Style.RESET_ALL}")
            print(f"  --> Fraction: {count/total_MC_events:.4f}\n")