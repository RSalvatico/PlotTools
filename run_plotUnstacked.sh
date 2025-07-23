#!/bin/sh

# Unstacked scores for signal, ttLF, total background
#python3 plotUnstacked.py --input_dir histos_07072025/scores_ttLFm0p1/ --input_csv hconfig_scores.csv --output_dir test/ --process ttLF --log
#python3 plotUnstacked.py --input_dir histos_07072025/scores_ttLFm0p1/ --input_csv hconfig_scores.csv --output_dir test/ --process ttbb --log
#python3 plotUnstacked.py --input_dir histos_07072025/scores_ttLFm0p1/ --input_csv hconfig_scores.csv --output_dir test/ --process ttbj --log
#python3 plotUnstacked.py --input_dir histos_07072025/scores_ttLFm0p1/ --input_csv hconfig_scores.csv --output_dir test/ --process ttcc --log
#python3 plotUnstacked.py --input_dir histos_07072025/scores_ttLFm0p1/ --input_csv hconfig_scores.csv --output_dir test/ --process ttcj --log

# Purity/evt number plots for CRs and SR
#python3 plotUnstacked.py --input_dir histos_03062025_scores/ --output_dir purity_plots/CRSR/ --purity --multiRegion
#python3 plotUnstacked.py --input_dir histos_03062025_scores/ --output_dir purity_plots/CRSR/ --purity --multiRegion --raw_evt_number

# Purity/evt number plots for 4FS vs 5FS, ttbb and ttbj
#python3 plotUnstacked.py --input_dir histos_07072025/ --output_dir purity_plots/FS/ --plot_4F5F --process ttbb
#python3 plotUnstacked.py --input_dir histos_07072025/ --output_dir purity_plots/FS/ --plot_4F5F --process ttbb --raw_evt_number
#python3 plotUnstacked.py --input_dir histos_07072025/ --output_dir purity_plots/FS/ --plot_4F5F --process ttbj
#python3 plotUnstacked.py --input_dir histos_07072025/ --output_dir purity_plots/FS/ --plot_4F5F --process ttbj --raw_evt_number

# 4FS vs 5FS comparison of ttbb+ttbj in the ttWcb score
python3 plotUnstacked.py --input_dir histos_07072025/ --output_dir purity_plots/FS_vs_score/ --plot_4F5F_vs_score