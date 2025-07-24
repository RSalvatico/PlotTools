import ROOT
import argparse
import glob
import math
import csv
import os
import cmsstyle as CMS

def stack_histograms(input_files, hist_name, output_dir, sonly, sig_norm, log, blind):
    """
    Reads TH1Ds with the same name from multiple files, stacks them in a THStack, and saves the result.

    Parameters:
    - input_dir: Input directory, where ROOT files are located.
    - hist_name: Name of the histograms to stack.
    - output_dir: Output directory for the TCanvas containing THStacks.
    - sonly: Decide whether to plot only the signal.
    - sig_norm: Normalization of the signal.
    - log: Use log scale on the Y-axis.
    - blind: Decide whether to blind the data in (b,c) invariant mass histogram.
    """
    # Create a THStack and a dictionary {process name : histogram} to feed to the CMS plotting
    stack = ROOT.THStack("stack", f"Stack of {hist_name}")
    phys_process = dict()

    # Define X-axis boundaries for the stack
    x_low, x_high = 0., 1.

    # Create a histogram for the signal and the collision data to be added separately from the stack
    sig_hist = ROOT.TH1D()
    data_hist = ROOT.TH1D()

    if sonly:
        print(f"Plotting only the signal")

    # Decide whether to blind the data in the invariant mass histogram
    isBlind = True if (hist_name == "h_mass_minDR_bc" and blind) else False

    # Open input files and retrieve histograms
    for infile in input_files:
        
        if sonly and "Wcb" not in infile:
            continue

        print(f"Reading file: {infile}")

        # Open the file
        root_file = ROOT.TFile.Open(infile)
        if not root_file or root_file.IsZombie():
            raise FileNotFoundError(f"Could not open file: {infile}")

        # Retrieve the histogram
        hist = root_file.Get(hist_name)
        if not hist or not isinstance(hist, ROOT.TH1):
            raise ValueError(f"Histogram '{hist_name}' not found in file '{infile}'.")

        # Clone the histogram to avoid issues when the file is closed
        hist_clone = hist.Clone()
        hist_clone.SetDirectory(0)  # Detach from the file

        # Assign X-axis boundaries for the stack
        x_low = hist_clone.GetBinLowEdge(1)
        x_high = hist_clone.GetBinLowEdge(hist_clone.GetNbinsX() + 1)

        # Treat the signal sample separately (it will be added also as a dashed line to the plots)
        if "Wcb" in infile:
            print(f"W->cb x {sig_norm} histogram will be also added to the plot separately")
            sig_hist = hist.Clone()
            sig_hist.SetDirectory(0)
            sig_hist *= sig_norm
            if sonly: continue # Avoid adding W->cb to the stack when plotting signal only
        if "Data" in infile:
            print(f"Data histogram will be added to the plot separately")
            if isBlind: continue
            data_hist = hist_clone
            continue

        # Fill dictionary {process name : histogram} to feed to the CMS plotting
        phys_process_name = (infile.split('_')[-1]).replace('.root','')
        phys_process[phys_process_name] = hist_clone

        # Close the file
        root_file.Close()

    # Save the stack in a canvas and add a legend
    print(f"Saving stacked histograms as: {output_dir}{hist_name.replace('h_','')}.pdf")
    canvas = CMS.cmsDiCanvas('canvas', x_low, x_high, 0, 1, 0.7, 1.3, hist_name.replace('h_',''), 'Events', 'Data/MC', square = CMS.kSquare, extraSpace=0.01, iPos=11)
    canvas.cd(1)
    legend = CMS.cmsLeg(0.65,0.4,0.85,0.87, textSize=0.04) # Needs to be defined after the cmsCanvas or it won't be plotted
    if not sonly and not isBlind:
        legend.AddEntry(data_hist, "Data", "pe")
    legend.AddEntry(sig_hist, f"W#rightarrow cb #times {sig_norm}", "l")

    CMS.cmsDrawStack(stack,legend,phys_process)
    if not sonly:
        CMS.cmsDraw(sig_hist,"same", lstyle = 2, msize = 0, lcolor = ROOT.kRed, lwidth = 4)
    else:
        CMS.cmsDraw(sig_hist,"same, hist", msize = 0, fcolor = ROOT.kRed, lcolor = ROOT.kRed, fstyle = 3018)
    CMS.cmsDraw(data_hist, "E1X0", mcolor=ROOT.kBlack)

    # Set Y-axis range based on maximum value of stacked histograms
    hist_from_canvas = CMS.GetcmsCanvasHist(canvas.GetPad(1))
    hist_from_canvas.GetYaxis().SetRangeUser(0.01,max(stack.GetHistogram().GetMaximum(),data_hist.GetMaximum()) * 1.2)
    if sonly:
        hist_from_canvas.GetYaxis().SetRangeUser(0.01,sig_hist.GetMaximum() * 1.2)
    hist_from_canvas.GetYaxis().SetMaxDigits(3) # Force scientific notation above 3 digits on the Y-axis
    #Draw the stack in log scale
    if log: 
        ROOT.gPad.SetLogy()
        hist_from_canvas.GetYaxis().SetRangeUser(0.0001,max(stack.GetHistogram().GetMaximum(),data_hist.GetMaximum()) * 10000)
        if sonly:
            hist_from_canvas.GetYaxis().SetRangeUser(0.01,sig_hist.GetMaximum() * 1000)

    # Add error bars
    if not sonly and not isBlind:
        bkg_hist = stack.GetStack().Last()
        err_hist = bkg_hist.Clone()
        CMS.cmsDraw(err_hist, "e2same0", lcolor = 335, lwidth = 1, msize = 0, fcolor = ROOT.kBlack, fstyle = 3004) 

        # Ratio plot
        canvas.cd(2)
        ratio = data_hist.Clone("ratio")
        ratio.Divide(bkg_hist)

        for i in range(1,ratio.GetNbinsX()+1):
            if(ratio.GetBinContent(i)):
                ratio.SetBinError(i, math.sqrt(data_hist.GetBinContent(i))/bkg_hist.GetBinContent(i))
            else:
                ratio.SetBinError(i, 10^(-99))

        yerr = ROOT.TGraphAsymmErrors()
        yerr.Divide(data_hist, bkg_hist, 'pois') 
        for i in range(0,yerr.GetN()+1):
            yerr.SetPointY(i,1)
        CMS.cmsDraw(yerr, "e2same0", lwidth = 100, msize = 0, fcolor = ROOT.kBlack, fstyle = 3004)  
        CMS.cmsDraw(ratio, "E1X0", mcolor=ROOT.kBlack)
        ref_line = ROOT.TLine(x_low, 1, x_high, 1)
        CMS.cmsDrawLine(ref_line, lcolor = ROOT.kBlack, lstyle = ROOT.kDotted)
        ratio_from_canvas = CMS.GetcmsCanvasHist(canvas.GetPad(2))
        ratio_from_canvas.GetYaxis().SetRangeUser(0.5,1.5)

    # Save the canvas in pdf and png formats
    plot_name = f"{output_dir}{hist_name.replace('h_','')}" if not log else f"{output_dir}/log/{hist_name.replace('h_','')}"
    CMS.SaveCanvas(canvas,f"{plot_name}.pdf") # The False is needed not to close the canvas
    print()

def create_output_dir(output_dir, log):
    """
    Create the output directory if it does not exist.

    Parameters:
    - output_dir: The directory where the output files will be saved.
    - log: Boolean indicating whether a log directory should be created.
    """
    os.makedirs(output_dir, exist_ok=True)
    if log:
        os.makedirs(os.path.join(output_dir, 'log'), exist_ok=True)

def read_csv(csv_file):
    """
    Open and read a csv file containing the name and the range of the variables to be plotted. 
    Fill in a list of histogram names.

    Parameters:
    - csv_file: The csv file containing variable names and binning for the respective histograms.
    """
    with open(csv_file, mode = 'r') as f:
        csv_reader = csv.reader(f) 
        hist_list = [line[0] for line in csv_reader if not line[0] == 'Variable']
        hist_list = [f"h_{hist}" for hist in hist_list]

    return hist_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stack TH1D histograms from multiple ROOT files.")
    parser.add_argument("--input_dir", type=str, required=True, help="Input directory, where ROOT files are located.")
    parser.add_argument("--hist_name", required=False, help="Name of the histograms to stack.")
    parser.add_argument("--input_csv", type=str, required=True, help="The csv file to read variables and ranges from.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for the TCanvas containing THStacks.")
    parser.add_argument("--sonly", nargs="?", const=1, type=bool, default=False, required=False, help="Decide whether to plot only the signal.")
    parser.add_argument("--sig_norm", nargs="?", const=1, type=int, default=1, required=False, help="Signal normalization.")
    parser.add_argument("--log", nargs="?", const=1, type=bool, default=False, required=False, help="Decide whether to use log scale on the Y-axis.")
    parser.add_argument("--blind", nargs="?", const=1, type=bool, default=False, required=False, help="Decide whether to blind the data in invarian mass histogram.")

    args = parser.parse_args()

    # Set plotting details
    CMS.SetExtraText("Work in progress")
    CMS.SetLumi("59.83")

    # Get input files from the input_dir
    input_files = glob.glob(f"{args.input_dir}*.root")
    # Sort the input files to ensure consistent ordering
    input_files.sort()

    # Create the output directory if it does not exist
    create_output_dir(args.output_dir, args.log)

    # Plot either all histograms from the csv file or a single histogram
    if not args.hist_name:
        hist_list = read_csv(args.input_csv)
        for hist_name in hist_list:
            stack_histograms(input_files, hist_name, args.output_dir, args.sonly, args.sig_norm, args.log, args.blind)
    else:
        stack_histograms(input_files, args.hist_name, args.output_dir, args.sonly, args.sig_norm, args.log, args.blind)
