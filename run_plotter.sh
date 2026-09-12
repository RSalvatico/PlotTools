#!/bin/sh
# Preselection score
PROD_VERSION=07042026
EXTRA_NAME=fscores_onlyEle
#EXTRA_NAME=preselection_withSysts/
INPUT_DIR=histos_$PROD_VERSION/$EXTRA_NAME/
OUTPUT_DIR=plots_$PROD_VERSION/$EXTRA_NAME/
EXTRA_PATH=fscores_onlyEle
#EXTRA_PATH=preselection_withSysts/
SIG_NORM=5
CONFIG_FILE=configs/hconfig_fscores.csv

python3 plotter.py --input_dir $INPUT_DIR --output_dir $OUTPUT_DIR --sig_norm $SIG_NORM --input_csv $CONFIG_FILE --blind
python3 plotter.py --input_dir $INPUT_DIR --output_dir $OUTPUT_DIR --sig_norm $SIG_NORM --input_csv $CONFIG_FILE --blind --log

#python3 plotter.py --input_dir histos_07012026/ttLFm0p1/allPlots/ --output_dir plots_07012026/ttLFm0p1_ttWcbm0p7/allPlots/ --sig_norm $SIG_NORM --input_csv configs/hconfig.csv --blind
#python3 plotter.py --input_dir histos_07012026/ttLFm0p1/allPlots/ --output_dir plots_07012026/ttLFm0p1_ttWcbm0p7/allPlots/ --sig_norm $SIG_NORM --input_csv configs/hconfig.csv --blind --log

#python3 plotter.py --input_dir histos_centralVcb_njets4/ --output_dir plots_centralVcb_njets4/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_njets5/ --output_dir plots_centralVcb_njets5/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_njets6/ --output_dir plots_centralVcb_njets6/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_njetsM6/ --output_dir plots_centralVcb_njetsM6/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nbtag3/ --output_dir plots_centralVcb_nbtag3/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nbtag4/ --output_dir plots_centralVcb_nbtag4/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nbtag5/ --output_dir plots_centralVcb_nbtag5/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nbtagM5/ --output_dir plots_centralVcb_nbtagM5/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nctag1/ --output_dir plots_centralVcb_nctag1/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nctag2/ --output_dir plots_centralVcb_nctag2/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nctag3/ --output_dir plots_centralVcb_nctag3/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log
#python3 plotter.py --input_dir histos_centralVcb_nctagM3/ --output_dir plots_centralVcb_nctagM3/ --sig_norm $SIG_NORM --input_csv configs/hconfig_scores.csv --log

# Scores
#python3 plotter.py --input_dir histos_centralVcb/fscores_ttLFm0p1_rebinned/ --output_dir plots_centralVcb/fscores_ttLFm0p1_rebinned/ --sig_norm $SIG_NORM --input_csv configs/hconfig_fscores.csv --blind
#python3 plotter.py --input_dir histos_centralVcb/fscores_ttLFm0p1_rebinned/ --output_dir plots_centralVcb/fscores_ttLFm0p1_rebinned/ --sig_norm $SIG_NORM --input_csv configs/hconfig_fscores.csv --blind --log

# Copy the plots to the appropriate web area
mkdir -p /eos/user/r/rselvati/www/Vcb/Run3/$PROD_VERSION/$EXTRA_PATH/log/
cp /eos/user/r/rselvati/www/index.php /eos/user/r/rselvati/www/Vcb/Run3/$PROD_VERSION/
cp /eos/user/r/rselvati/www/index.php /eos/user/r/rselvati/www/Vcb/Run3/$PROD_VERSION/$EXTRA_PATH/
cp /eos/user/r/rselvati/www/index.php /eos/user/r/rselvati/www/Vcb/Run3/$PROD_VERSION/$EXTRA_PATH/log/
cp $OUTPUT_DIR/* /eos/user/r/rselvati/www/Vcb/Run3/$PROD_VERSION/$EXTRA_PATH/
cp $OUTPUT_DIR/log/* /eos/user/r/rselvati/www/Vcb/Run3/$PROD_VERSION/$EXTRA_PATH/log/