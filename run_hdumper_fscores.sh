#!/bin/sh
# This script is used to run the hdumper to make histograms with fscores (i.e., those to fit)
python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/mc/ --output_dir histos_07072025/4F/  --tree_name Events --input_csv hconfig_fscores.csv --year 2018 --eventClassification
python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/data/ --output_dir histos_07072025/4F/  --tree_name Events --input_csv hconfig_fscores.csv --year 2018 --eventClassification

#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/mc/ --output_dir histos_07072025/  --tree_name Events --input_csv hconfig_fscores.csv --year 2018 --eventClassification --use5FS
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/data/ --output_dir histos_07072025/  --tree_name Events --input_csv hconfig_fscores.csv --year 2018 --eventClassification --use5FS