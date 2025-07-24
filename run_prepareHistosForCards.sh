#!/bin/sh
#python3 prepareHistosForCards.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/data/total/ --output_dir test_withSyst_newHists/ --tree_name Events --year 2018
#python3 prepareHistosForCards.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/mc/ --output_dir test_withSyst_newHists/ --tree_name Events --year 2018

python3 prepareHistosForCards.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/data/total/ --output_dir datacard_preparation/ --tree_name Events --year 2018
python3 prepareHistosForCards.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/29012025_2018_1L/mc/ --output_dir datacard_preparation/ --tree_name Events --year 2018