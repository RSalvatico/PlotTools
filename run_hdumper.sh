#!/bin/sh
# This script is used to run the hdumper to make histograms
python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/mc/ --output_dir histos_07072025/scores_ttLFm0p1_5FS/  --tree_name Events --input_csv hconfig_scores.csv --year 2018 --add_selection "score_ttLF < 0.1" --use5FS
python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/data/ --output_dir histos_07072025/scores_ttLFm0p1_5FS/  --tree_name Events --input_csv hconfig_scores.csv --year 2018 --add_selection "score_ttLF < 0.1" --use5FS