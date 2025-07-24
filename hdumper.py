import ROOT
import argparse
import glob
import csv
import os
import sys
from colorama import Fore, Style 
import numpy as np

ROOT.ROOT.EnableImplicitMT()

def process_trees(input_files, output_files, tree_name, hist_configs, year, selections, eventClassification, use5FS):
    """
    Processes multiple TTrees, converts them to multiple TH1Ds for specified branches, and saves them to ROOT files.

    Parameters:
    - input_files: List of input ROOT files.
    - output_files: List of output ROOT files.
    - tree_names: List of TTree names corresponding to input files.
    - hist_configs: List of dictionaries with keys 'branch', 'nbins', 'xmin', 'xmax'.
    - year: Year of data taking.
    - selections: String containing common event preselection.
    - eventClassification: Boolean indicating whether to apply event classification.
    - use5FS: Boolean indicating whether to use 5-flavor scheme MC for ttbb and ttbj processes.
    """
    print("")
    if not (len(input_files) == len(output_files)):
        raise ValueError("Input files and output files must have the same length.")

    for infile, outfile in zip(input_files, output_files):

        print(f"{Fore.RED}Processing file: {infile}{Style.RESET_ALL}")

        # Open input file
        input_file = ROOT.TFile.Open(infile)
        if not input_file or input_file.IsZombie():
            raise FileNotFoundError(f"Could not open file: {infile}")

        # Access the TTree
        tree = input_file.Get(tree_name)
        if not tree or not isinstance(tree, ROOT.TTree):
            raise ValueError(f"TTree '{tree_name}' not found in file '{infile}'.")

        # Create RDataFrame from TTree
        df = ROOT.RDataFrame(tree)

        if eventClassification:
            print(f"{Fore.YELLOW}Running in event classification mode. Will define a series of fractional scores.{Style.RESET_ALL}")
            # Define the fractional scores
            df = df.Define("fscore_ttbb", "score_ttbb / (score_ttbb + score_ttbj + score_ttcc + score_ttcj + score_ttLF)") \
                .Define("fscore_ttbj", "score_ttbj / (score_ttbb + score_ttbj + score_ttcc + score_ttcj + score_ttLF)") \
                .Define("fscore_ttcc", "score_ttcc / (score_ttbb + score_ttbj + score_ttcc + score_ttcj + score_ttLF)") \
                .Define("fscore_ttcj", "score_ttcj / (score_ttbb + score_ttbj + score_ttcc + score_ttcj + score_ttLF)") \
                .Define("fscore_ttLF", "score_ttLF / (score_ttbb + score_ttbj + score_ttcc + score_ttcj + score_ttLF)")
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

        tt_file_names = ["ttbb-4f", "ttbb-dps", "ttbar-powheg"]
        tt4f_strings = ["ttbb", "ttbj"]
        tt_strings   = ["ttcc", "ttcj", "ttLF"]

        # Assign event weight based on data taking year and process type
        weight = assign_event_weight(year, infile)

        # Process each selection-output combinations
        for selection_name in selections:

            # Apply base selection to every sample; apply the ttbar-specific selection to the right 4F, dps, and 5F powheg samples
            if not "base" in selection_name and not any(x in infile for x in tt_file_names): 
                continue
            if any(x in infile for x in tt_file_names) and "base" in selection_name:
                continue
            if use5FS: # "ttbb", "ttbj" -> both powheg and dps samples; "ttcc", "ttcj", "ttLF" --> only powheg
                if any(x in selection_name for x in tt4f_strings) and not ("powheg" in infile or "dps" in infile):
                    continue
                if any(x in selection_name for x in tt_strings) and not "powheg" in infile: 
                    continue
            else:
                if any(x in selection_name for x in tt4f_strings) and not "bb" in infile:
                    continue
                if any(x in selection_name for x in tt_strings) and not "powheg" in infile:
                    continue

        
            # Produce output file
            tt_outfile_name = outfile.replace('.root','_'+selection_name+'.root')
            output_file = tt_outfile_name if not "base" in selection_name else outfile
            output_root = ROOT.TFile(output_file, "RECREATE")

            # Add event selection making sure that the "base" selection is applied everywhere
            event_selection = f"{selections['base']}{selections[selection_name]}" if not "base" in selection_name else f"{selections[selection_name]}"
            if "singlee" in infile:
                event_selection += " && passTrigMu==0" # Remove from the electron channel the events that fired the muon trigger. Could choose to do vice versa as well.
            print(f"Applying selection: {Fore.GREEN}{event_selection}{Style.RESET_ALL} -> Producing output file: {output_file}")
            df_selected = df.Filter(event_selection)

            #print("after first event selection:", df_selected.Count().GetValue())

            # If weight is a complex expression, define it as a new column
            weight_column = "weight_column"
            if not "data" in infile:
                print(f"Event weight: {weight}")
                df_selected = df_selected.Define(weight_column, weight)
            else: # Keep the weight 1 for collision data
                df_selected = df_selected.Define(weight_column, "1")

            # Define event classification for the dedicated mode
            if eventClassification:
                from weights_and_constants import adhoc_selection, adhoc_binning
                adhoc_selection = adhoc_selection.copy()
                adhoc_binning = adhoc_binning.copy()
                #eventClassificationBaseSelection = "score_tt_Wcb > 0.6 && score_ttLF < 0.05"
                #SR_selection = "score_tt_Wcb > 0.85"
                #CR_selection = "score_tt_Wcb < 0.85"
                #adhoc_selection = {
                #    "score_tt_Wcb" : f"{eventClassificationBaseSelection} && {SR_selection}",
                #    "fscore_ttbb"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttcj']} * score_ttcj && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttLF']} * score_ttLF",
                #    "fscore_ttbj"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttcj']} * score_ttcj && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttLF']} * score_ttLF",
                #    "fscore_ttcc"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttcj']} * score_ttcj && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttLF']} * score_ttLF",
                #    "fscore_ttcj"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttLF']} * score_ttLF",
                #    "fscore_ttLF"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttcj']} * score_ttcj"
                #}
                #adhoc_binning = {
                #    "score_tt_Wcb" : np.array([0.,0.9,1.]),
                #    "fscore_ttbb"  : np.array([0.,0,0.7,1.]),
                #    "fscore_ttbj"  : np.array([0.,0.45,1.]),
                #    "fscore_ttcc"  : np.array([0.,0.45,1.]),
                #    "fscore_ttcj"  : np.array([0.,0.35,1.]),
                #    "fscore_ttLF"  : np.array([0.,0.1,1.]),
                #}

            # Create histograms for each branch
            final_df = dict()
            for hist_config in hist_configs:
                branch_name = hist_config['branch']
                nbins = int(hist_config['nbins'])
                xmin = float(hist_config['xmin'])
                xmax = float(hist_config['xmax'])
                print(f"Creating histogram for branch: {branch_name}")
                
                final_df[branch_name] = df_selected.Filter(adhoc_selection[branch_name]) if eventClassification else df_selected

                # Create histogram
                if eventClassification:
                    hist = final_df[branch_name].Histo1D((f"h_{branch_name}", f"Histogram of {branch_name}", len(adhoc_binning[branch_name])-1, adhoc_binning[branch_name]), branch_name, weight_column)
                else:
                    hist = final_df[branch_name].Histo1D((f"h_{branch_name}", f"Histogram of {branch_name}", nbins, xmin, xmax), branch_name, weight_column)

                #print(f"Number of events after full selection: {hist.Integral()}")


                # Write histogram to output file
                hist.Write()

            # Close files
            output_root.Close()

        input_file.Close()

        print(f"Saved histograms to: {outfile}\n")

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

def assign_event_weight(year, infile):
    """
    Define the MC event weight according to the year. Collision data should be handled separately.

    Parameters:
    - year: Data taking year.
    - infile: Input file.
    """
    weight = "1"
    if year == 2018:
        weight = "lumiwgt*genWeight*xsecWeight*l1PreFiringWeight*puWeight*muEffWeight*elEffWeight*flavTagWeight*(((abs(lep1_pdgId)==11 && passTrigEl && ((year!=2018) || (year==2018 && !(lep1_phi>-1.57 && lep1_phi<-0.87 && lep1_eta<-1.3)))) || (abs(lep1_pdgId)==13 && passTrigMu)) && passmetfilters)"
    if "ttbar" in infile:
        weight = f"{weight}*topptWeight"
    if "4f" in infile:
        weight = f"{weight}*topptWeight"#*0.7559" # 5FS / 4FS for tt+B component
    
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

def score_calculation(score_tt_Wcb, score_ttLF, score_ttbb, score_ttbj):
    """
    Calculate the fractional score for the ttbb, ttbj, ttcc, and ttLF samples.

    Parameters:
    - scores corresponding to the different processes.
    """
    return score_tt_Wcb / (score_tt_Wcb + score_ttbb + score_ttbj + score_ttLF)


    
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
    parser.add_argument("--event_counting_file", type=str, required=False, help="File to save event counts for each selection.")
    parser.add_argument("--eventClassification", nargs="?", const=1, type=bool, default=False, required=False, help="Apply event classification selection.")
    parser.add_argument("--use5FS", nargs="?", const=1, type=bool, default=False, required=False, help="Use 5-flavor scheme.")

    args = parser.parse_args()

    # Get input files from the input_dirs list
    for input_dir in args.input_dirs:
        input_files = glob.glob(f"{input_dir}*.root")

    # Prepare list of output files based on the name of the input files
    output_files = prepare_output(args.output_dir, input_files)

    # Prepare histogram configurations for each branch
    hist_configs = read_csv(args.input_csv)

    selections = {"base": "n_ak4>=4 && (n_btagM+n_ctagM)>=3 && n_btagM>=1",
                 "ttbb" : " && genEventClassifier==9 && wcb==0",
                 "ttbj" : " && (genEventClassifier==7 || genEventClassifier==8) && wcb==0",
                 "ttcc" : " && genEventClassifier==6 && wcb==0",
                 "ttcj" : " && (genEventClassifier==4 || genEventClassifier==5) && wcb==0",
                 "ttLF" : " && tt_category==0 && higgs_decay==0 && wcb==0"
    }

    # These weights correspond (roughly) to the fraction of events of a certain process expected in the corresponding category after the ttWcb and ttLF score selection.
    #from weights_and_constants import weights_and_constants
    #wc = weights_and_constants()
    #evtClassification_weights = wc.weights_0p6ttWcb_and_0p1ttLF

    # Apply trigger selection to separate channels if requested
    if args.electron:
        selections["base"] += " && passTrigEl"
    if args.muon:
        selections["base"] += " && passTrigMu"

    use5FS = False
    if args.use5FS:
        use5FS = True
        print(f"{Fore.GREEN}Using 5-flavor scheme for ttbb and ttbj processes.{Style.RESET_ALL}")

    # Apply additional selections if specified
    if args.add_selection:
        for key in selections.keys():
            selections[key] += f" && ({args.add_selection})"

    process_trees(input_files, output_files, args.tree_name, hist_configs, args.year, selections, args.eventClassification, use5FS)

    # Merge some of the output files
    ttV_list = ["h_ttW.root", "h_ttZ.root"]
    merge_files(args.output_dir, ttV_list, "h_ttV.root")
    ttH_list = ["h_ttHbb.root", "h_ttHcc.root", "h_ttV.root"]
    merge_files(args.output_dir, ttH_list, "h_ttH-ttV.root")
    if use5FS:
        ttbb_list = ["h_ttbar-powheg_ttbb.root", "h_ttbb-dps_ttbb.root"]
        merge_files(args.output_dir, ttbb_list, "h_ttbb-withDPS.root")
        ttbj_list = ["h_ttbar-powheg_ttbj.root", "h_ttbb-dps_ttbj.root"]
        merge_files(args.output_dir, ttbj_list, "h_ttbj-withDPS.root")
    else:
        ttbb_list = ["h_ttbb-4f_ttbb.root", "h_ttbb-dps_ttbb.root"]
        merge_files(args.output_dir, ttbb_list, "h_ttbb-withDPS.root")
        ttbj_list = ["h_ttbb-4f_ttbj.root", "h_ttbb-dps_ttbj.root"]
        merge_files(args.output_dir, ttbj_list, "h_ttbj-withDPS.root")
    diboson_list = ["h_TWZ.root", "h_diboson.root"]
    merge_files(args.output_dir, diboson_list, "h_diboson-tWZ.root")
    data_list = ["h_singlee.root", "h_singlemu.root"]
    merge_files(args.output_dir, data_list, "h_Data.root")
