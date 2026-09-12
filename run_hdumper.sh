#!/bin/sh
# This script is used to run the hdumper to make histograms
INPUT_DIR=/eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07042026_2024_1L_Wcb/
PROD_VERSION=07042026
CONFIG_FILE=configs/hconfig_minimal.csv
#EXTRA_NAME=preselection_ge2bge1c_new2024FlavTagSFs

YEAR=2024

declare -A w=(
    [ttbb]=0.0709 [tt2b]=0.0378 [ttbj]=0.0863
    [ttcc]=0.0676 [tt2c]=0.0868 [ttcj]=0.1154 [ttLF]=0.5351
)
categories=(ttbb tt2b ttbj ttcc tt2c ttcj ttLF)

cat=ttbj
EXTRA_NAME=CR_${cat}_selection_ge2bge1c_lepEta_rescaledDPS_corrected
OUTPUT_DIR=histos_$PROD_VERSION/$EXTRA_NAME/
parts=()
for other in "${categories[@]}"; do
    [[ "$other" == "$cat" ]] && continue
    parts+=("${w[$cat]}*score_${cat} > ${w[$other]}*score_${other}")
done
conditions=$(printf '%s && ' "${parts[@]}"); conditions=${conditions% && }

EXTRA_SELECTION="score_ttLF < 0.1 && score_tt_Wcb < 0.8"
OUTPUT_DIR=histos_$PROD_VERSION/$EXTRA_NAME/
FLAVTAG_SF_JSON=$CMSSW_BASE/src/PhysicsTools/NanoTTH/data/flavTagSF/flavTaggingSF_2024.json.gz

python3 hdumper.py --input_dirs $INPUT_DIR/mc/ --output_dir $OUTPUT_DIR --tree_name Events --input_csv $CONFIG_FILE --year $YEAR --add_selection "$EXTRA_SELECTION && $conditions" --flavtag_sf_json $FLAVTAG_SF_JSON
python3 hdumper.py --input_dirs $INPUT_DIR/data/ --output_dir $OUTPUT_DIR --tree_name Events --input_csv $CONFIG_FILE --year $YEAR --add_selection "$EXTRA_SELECTION && $conditions" --flavtag_sf_json $FLAVTAG_SF_JSON

#python3 hdumper.py --input_dirs $INPUT_DIR/mc/ --output_dir $OUTPUT_DIR --tree_name Events --input_csv $CONFIG_FILE --year $YEAR --add_selection "$EXTRA_SELECTION && $conditions" --systematics --use5FS
#python3 hdumper.py --input_dirs $INPUT_DIR/data/ --output_dir $OUTPUT_DIR --tree_name Events --input_csv $CONFIG_FILE --year $YEAR --add_selection "$EXTRA_SELECTION && $conditions" --systematics --use5FS

#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1_ttWcbm0p7/allPlots/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "score_ttLF < 0.1 && score_tt_Wcb < 0.7"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/data/ --output_dir histos_07012026/ttLFm0p1_ttWcbm0p7/allPlots/ --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "score_ttLF < 0.1 && score_tt_Wcb < 0.7"

# Select ttLF < 0.1
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_18122025_promptMVA_2024_1L_Wcb/mc/ --output_dir histos_18122025_promptMVA/ttLFm0p1/score_tt_Wcb/4FS/  --tree_name Events --input_csv configs/hconfig_minimal.csv --year 2024 --add_selection "score_ttLF < 0.1"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_18122025_promptMVA_2024_1L_Wcb/mc/ --output_dir histos_18122025_promptMVA/ttLFm0p1/score_tt_Wcb/5FS/  --tree_name Events --input_csv configs/hconfig_minimal.csv --year 2024 --add_selection "score_ttLF < 0.1" --use5FS

# Plotting of the scores in the CR (no fscores defined) or in SR
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/scores/  --tree_name Events --input_csv configs/hconfig_scores.csv --year 2024
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/scores/ttLFm0p1/  --tree_name Events --input_csv configs/hconfig_scores.csv --year 2024 --add_selection "score_ttLF < 0.1"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1/CR/  --tree_name Events --input_csv configs/hconfig_scores.csv --year 2024 --add_selection "score_ttLF < 0.1 && score_tt_Wcb < 0.8"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1/SR/  --tree_name Events --input_csv configs/hconfig_scores.csv --year 2024 --add_selection "score_ttLF < 0.1 && score_tt_Wcb > 0.8"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1/CRfscores/  --tree_name Events --input_csv configs/hconfig_fscores.csv --year 2024 --eventClassification

# Plotting of the 4FS and 5FS (f)scores
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1/4FS/  --tree_name Events --input_csv configs/hconfig_fscores.csv --year 2024 --eventClassification
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1/5FS/  --tree_name Events --input_csv configs/hconfig_fscores.csv --year 2024 --eventClassification --use5FS
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1/score_tt_Wcb/4FS/  --tree_name Events --input_csv configs/hconfig_minimal.csv --year 2024 --add_selection "score_ttLF < 0.1"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_07012026_2024_1L_Wcb/mc/ --output_dir histos_07012026/ttLFm0p1/score_tt_Wcb/5FS/  --tree_name Events --input_csv configs/hconfig_minimal.csv --year 2024 --add_selection "score_ttLF < 0.1" --use5FS

# Plot a bit of everything in different n_ak4, n_btag and n_ctag bins
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_njets4/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 == 4"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_njets4/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 == 4"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_njets5/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 == 5"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_njets5/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 == 5"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_njets6/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 == 6"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_njets6/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 == 6"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_njetsM6/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 > 6"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_njetsM6/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ak4 > 6"

#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nbtag3/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM == 3"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nbtag3/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM == 3"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nbtag4/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM == 4"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nbtag4/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM == 4"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nbtag5/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM == 5"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nbtag5/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM == 5"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nbtagM5/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM > 5"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nbtagM5/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_btagM > 5"

#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nctag1/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM == 1"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nctag1/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM == 1"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nctag2/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM == 2"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nctag2/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM == 2"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nctag3/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM == 3"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nctag3/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM == 3"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/mc/ --output_dir histos_centralVcb_nctagM3/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM > 3"
#python3 hdumper.py --input_dirs /eos/cms/store/cmst3/group/top/rsalvatico/Vcb_analysis_centralVcb_2024_1L_Wcb/data/ --output_dir histos_centralVcb_nctagM3/  --tree_name Events --input_csv configs/hconfig.csv --year 2024 --add_selection "n_ctagM > 3"
