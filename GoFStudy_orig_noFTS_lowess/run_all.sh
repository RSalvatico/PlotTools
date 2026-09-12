#!/bin/bash
# Single entry point for the orig-shapes / lowess / no-flavtag-mirror study:
# CR and SR, expected (Asimov) and observed, prepostfit + impacts for all 4,
# plus goodness-of-fit toys for the two observed configs (CR_observed,
# SR_observed -- the expected/Asimov configs' toys were an earlier
# cross-check and are not part of this deliverable).
#
#   cmsenv && source setup.sh && source setup_rabbit.sh
#   cd GoFStudy_orig_noFTS_lowess && ./run_all.sh
#
# Idempotent: every step checks its own output file(s) and skips if already
# there. Safe to Ctrl-C and rerun, safe to rerun after a crash (see the AFS
# quota post-mortem in run_gof_toys_lean.sh's header) -- nothing here holds
# more than one toy-batch's hdf5 (~360 MB) on disk at a time.
#
# Toy count/parallelism: JOBS=16 TOYS=5000 ./run_all.sh to override; each
# rabbit fit is single-core, so JOBS should not exceed the machine's cores.
# The CR/SR toy batches alone take ~2-2.5h each at 16-way parallelism.
#
# BLINDING:
# - SR_observed unblinds tt-vcb IN THE FIT (real answer, real postfit shapes,
#   real GoF) -- the discipline is on reporting, not computation.
# - Its impacts/pulls plots pass --redact PARAM=tt-vcb to
#   rabbitPlotImpacts.py for every POI, not only when tt-vcb is the POI
#   being plotted: tt-vcb is a free parameter correlated with the tt+X norms,
#   so its real pull can appear as a text row on ANOTHER poi's impact plot
#   too (found and fixed once already -- see git log for the regression).
# - prepostfit (rabbitPlotStack.py) never draws an SR data marker unless
#   --show-signal-data is passed; SR_observed's default plots stay blinded,
#   a second _withSRdata copy is made for anyone who has cleared to see it.
# - This script never prints tt-vcb's fitted value to the terminal.
#
# GoF toy conditioning (see COMBINE_FIXES.md / run_gof_toys_lean.sh header):
# Asimov configs toy around the model's own prior (already the correct null
# for an Asimov check); observed configs toy around a postfit-conditioned
# copy of the tensor (built below), or Poisson noise gets double counted.

set -e

BASE=$(cd "$(dirname "$0")" && pwd)
cd "${BASE}"
#O=rabbit_5FS
O=rabbit
mkdir -p ${O} plots/prepostfit plots/impacts

#SHAPES=../datacards_07042026_correctTTweights_5FS/Vcb_SL_2024_shapes.root
SHAPES=../datacards_07042026_correctTTweights/Vcb_SL_2024_shapes.root
PM="--paramModel Mu --paramModel analysis.rabbit_models.FreeNorm ttbb,ttbj,tt2b,ttcc,ttcj,tt2c,ttLF"
FULL="--doImpacts --globalImpacts --saveHists --saveHistsPerProcess --computeHistErrors"
POIS="tt-vcb xsec_ttbb xsec_ttbj xsec_tt2b xsec_ttcc xsec_ttcj xsec_tt2c xsec_ttLF"
JOBS=${JOBS:-16}
TOYS=${TOYS:-5000}
TOYS_PER_JOB=$(( (TOYS + JOBS - 1) / JOBS ))
NPARAMS_IMPACTS=1000

# ------------------------------------------------------------------ [1/7] tensors
echo "### [1/7] tensors (lowess smoothing + flavTag mirror-up off -- prepareTensor.py defaults)"
[ -f ${O}/ourSR.hdf5 ] || python3 ../analysis/prepareTensor.py ${SHAPES} \
    -o ${O} --outname ourSR > ${O}/log_tensor_ourSR.log 2>&1 &
[ -f ${O}/ourCR.hdf5 ] || python3 ../analysis/prepareTensor.py ${SHAPES} \
    -o ${O} --outname ourCR --mask '_SR$' > ${O}/log_tensor_ourCR.log 2>&1 &
wait
echo "  done"

# ------------------------------------------------------------------ [2/7] fits
echo "### [2/7] fits: CR_expected, CR_observed, SR_expected (+tt-vcb scan), SR_observed"
[ -f ${O}/CR_expected.hdf5 ] || rabbit_fit.py ${O}/ourCR.hdf5 -o ${O} --outname CR_expected \
    -t -1 ${PM} --freezeParameters tt-vcb ${FULL} > ${O}/log_CR_expected.log 2>&1
[ -f ${O}/CR_observed.hdf5 ] || rabbit_fit.py ${O}/ourCR.hdf5 -o ${O} --outname CR_observed \
    -t 0 ${PM} --freezeParameters tt-vcb --unblind 'xsec_.*' ${FULL} > ${O}/log_CR_observed.log 2>&1
[ -f ${O}/SR_expected.hdf5 ] || rabbit_fit.py ${O}/ourSR.hdf5 -o ${O} --outname SR_expected \
    -t -1 ${PM} ${FULL} --scan tt-vcb --scanPoints 31 > ${O}/log_SR_expected.log 2>&1
if [ ! -f ${O}/SR_observed.hdf5 ]; then
    rabbit_fit.py ${O}/ourSR.hdf5 -o ${O} --outname SR_observed \
        -t 0 ${PM} --unblind 'tt-vcb' 'xsec_.*' ${FULL} > ${O}/log_SR_observed.log 2>&1
    grep -q "Unblinding 8 parameters:" ${O}/log_SR_observed.log \
        || { echo "ABORT: expected 8 unblinded parameters (tt-vcb + 7 xsec_*)" >&2; exit 1; }
fi
echo "  done"

# ------------------------------------------------------------------ [3/7] prepostfit
echo "### [3/7] prepostfit (5 = 4 configs + SR_observed's SR-data-visible twin)"
[ -f plots/prepostfit/CR_expected_postfit.png ] || python3 ../analysis/rabbitPlotStack.py \
    ${O}/CR_expected.hdf5 -o plots/prepostfit --postfix CR_expected --logy --asimov > /dev/null 2>&1
[ -f plots/prepostfit/CR_observed_postfit.png ] || python3 ../analysis/rabbitPlotStack.py \
    ${O}/CR_observed.hdf5 -o plots/prepostfit --postfix CR_observed --logy > /dev/null 2>&1
[ -f plots/prepostfit/SR_expected_postfit.png ] || python3 ../analysis/rabbitPlotStack.py \
    ${O}/SR_expected.hdf5 -o plots/prepostfit --postfix SR_expected --logy --asimov > /dev/null 2>&1
[ -f plots/prepostfit/SR_observed_postfit.png ] || python3 ../analysis/rabbitPlotStack.py \
    ${O}/SR_observed.hdf5 -o plots/prepostfit --postfix SR_observed --logy > /dev/null 2>&1
[ -f plots/prepostfit/SR_observed_withSRdata_postfit.png ] || python3 ../analysis/rabbitPlotStack.py \
    ${O}/SR_observed.hdf5 -o plots/prepostfit --postfix SR_observed_withSRdata --logy --show-signal-data \
    > /dev/null 2>&1
echo "  done"

# ------------------------------------------------------------------- [4/7] impacts
echo "### [4/7] impacts (traditional + global, 8 POIs x 4 configs = 64 plots)"
for CFG in CR_expected CR_observed SR_expected SR_observed; do
    [ -f plots/impacts/impacts_global_ungrouped_xsec_ttLF_${CFG}.png ] && continue
    ASI=()
    { [ "${CFG}" = "CR_expected" ] || [ "${CFG}" = "SR_expected" ]; } && ASI=(--asimov)
    REDACT=()
    [ "${CFG}" = "SR_observed" ] && REDACT=(--redact tt-vcb)
    for POI in ${POIS}; do
        for IT in traditional global; do
            python3 ../analysis/rabbitPlotImpacts.py ${O}/${CFG}.hdf5 --poi ${POI} \
                -o plots/impacts --postfix ${CFG} --impact-type ${IT} -n ${NPARAMS_IMPACTS}\
                "${ASI[@]}" "${REDACT[@]}" > /dev/null 2>&1
        done
    done
done
echo "  done"

# ------------------------------------------------------- [5/7] postfit-mean tensors
# CR_observed and SR_observed's own postfit expectation, patched into a copy
# of the tensor's data_obs, so --toysDataMode observed doesn't double count
# the real Poisson noise already baked into the raw data.
echo "### [5/7] postfit-conditioned tensors for GoF toys"
patch_postfit_mean() {  # <src tensor> <fitresult> <dst tensor>
    local SRC=$1 FIT=$2 DST=$3
    [ -f "${DST}" ] && return 0
    python3 - "${SRC}" "${FIT}" "${DST}" <<'PY'
import shutil, sys, h5py
import numpy as np
from rabbit import io_tools, inputdata

src, fit, dst = sys.argv[1], sys.argv[2], sys.argv[3]
indata = inputdata.FitInputData(src)
raw = np.asarray(indata.data_obs)

fr, meta = io_tools.get_fitresult(fit, meta=True)
chans = fr['mappings']['BaseMapping']['channels']

check = np.zeros_like(raw)
postfit = np.zeros_like(raw)
for ch, info in indata.channel_info.items():
    # Masked channels (e.g. the SR in the CR-only fit) carry no data_obs in
    # either the tensor's data_obs vector or the fitresult's per-channel
    # hists -- they are pure prediction. Skip them; there is nothing to
    # condition and chans[ch] has no 'hist_data_obs' key.
    if info.get('masked', False):
        continue
    s, e = info['start'], info['stop']
    check[s:e] = chans[ch]['hist_data_obs'].get().values().flatten()
    postfit[s:e] = chans[ch]['hist_postfit_inclusive'].get().values().flatten()
assert np.allclose(raw, check, rtol=0, atol=1e-6), "channel ordering mismatch, refusing to patch"
assert np.all(postfit > 0), "non-positive postfit bin, unsafe for Poisson toys"

shutil.copyfile(src, dst)
with h5py.File(dst, 'r+') as f:
    f['hdata_obs'][...] = postfit.astype(f['hdata_obs'].dtype)
print(f'wrote {dst}: sums raw={raw.sum():.1f} postfit={postfit.sum():.1f}')
PY
}
patch_postfit_mean ${O}/ourCR.hdf5 ${O}/CR_observed.hdf5 ${O}/ourCR_postfitmean.hdf5
patch_postfit_mean ${O}/ourSR.hdf5 ${O}/SR_observed.hdf5 ${O}/ourSR_postfitmean.hdf5
echo "  done"

# ------------------------------------------------------------------- [6/7] GoF toys
echo "### [6/7] GoF toys, CR_observed + SR_observed only (${TOYS} each, storage-safe batches)"
[ -f ${O}/toyGoF_CR_observed_summary.npz ] || JOBS=${JOBS} TOYS=${TOYS} SEED_BASE=60000 \
    ./run_gof_toys_lean.sh CR_observed ourCR_postfitmean.hdf5 observed \
    "--freezeParameters tt-vcb --unblind xsec_.*"
[ -f ${O}/toyGoF_SR_observed_summary.npz ] || JOBS=${JOBS} TOYS=${TOYS} SEED_BASE=50000 \
    ./run_gof_toys_lean.sh SR_observed ourSR_postfitmean.hdf5 observed \
    "--unblind tt-vcb xsec_.*"

# --------------------------------------------------------------------- [7/7] plots
echo "### [7/7] toy GoF plots"
for CFG in CR_observed SR_observed; do
    [ -f ${O}/toyGoF_${CFG}_summary.npz ] || continue
    python3 ../analysis/plotToyGoF.py --npz ${O}/toyGoF_${CFG}_summary.npz \
        --observed ${O}/${CFG}.hdf5 -o plots/toyGoF_${CFG}.png --label "Lepton + jets" \
        > /dev/null 2>&1
    echo "  plots/toyGoF_${CFG}.png"
done

echo "### ALL DONE -> plots/prepostfit, plots/impacts, plots/toyGoF_*.png"
