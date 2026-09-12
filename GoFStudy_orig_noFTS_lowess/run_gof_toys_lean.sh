#!/bin/bash
# GoF toy batches, storage-safe version.
#
# Post-mortem: run_bulk.sh's run_toy_batch() kept the full per-toy fit-result
# hdf5 (covariance, parms histogram, everything) for a study that only ever
# reads two scalars per toy (nllvalreduced, ndfsat). Measured cost ~1.1 MB/toy
# -- a 313-toy batch is ~360 MB, 16 batches/config ~5.7 GB. Four configs would
# need ~23 GB; the AFS volume quota (100 GB, work.sewuchte) hit 100% partway
# through CR_observed, corrupting all 16 of its batches (6 from the write
# failure itself, 10 more from having to kill them mid-write to stop the
# bleeding) and blocking SR_expected before a single toy ran.
#
# Fix: EACH batch, right after its own -t N run finishes, gets compacted by
# compact_toys.py (extract q=2*nllvalreduced + ndf into a tiny .npz, verify
# the round-trip, THEN delete the source .hdf5) before the next batch's slot
# opens up. Steady-state disk use for N toys/job is one job's hdf5 (~360 MB)
# at a time, not all 16 at once (~5.7 GB) -- an ~16x reduction in the peak.
#
#   cmsenv && source setup.sh && source setup_rabbit.sh
#   cd GoFStudy_orig_noFTS_lowess
#   ./run_gof_toys_lean.sh CR_observed ourCR_postfitmean.hdf5 observed \
#       "--freezeParameters tt-vcb --unblind xsec_.*"
#   ./run_gof_toys_lean.sh SR_expected ourSR.hdf5 expected ""
#
# TENSOR is resolved relative to the output dir O (default "rabbit"; override
# with e.g. O=rabbit_5FS ./run_gof_toys_lean.sh ... to run against a
# different tensor/fit set, such as the 5FS-only ttbar datacard).

set -e

NAME=$1 TENSOR=$2 MODE=$3 EXTRA=$4
[ -z "${NAME}" ] && { echo "usage: $0 <name> <tensor> <expected|observed> <extra fit flags>" >&2; exit 1; }

BASE=$(cd "$(dirname "$0")" && pwd)
cd "${BASE}"
O=${O:-rabbit}
PM="--paramModel Mu --paramModel analysis.rabbit_models.FreeNorm ttbb,ttbj,tt2b,ttcc,ttcj,tt2c,ttLF"
JOBS=${JOBS:-16}
TOYS=${TOYS:-5000}
TOYS_PER_JOB=$(( (TOYS + JOBS - 1) / JOBS ))
SEED_BASE=${SEED_BASE:-60000}

echo "### ${NAME}: ${JOBS} jobs x ${TOYS_PER_JOB} toys, compact-and-delete per batch"
for k in $(seq 1 ${JOBS}); do
    (
        N=toyGoF_${NAME}_batch${k}
        SUM=${O}/${N}_summary.npz
        [ -f "${SUM}" ] && exit 0   # already compacted, nothing to do
        rabbit_fit.py ${O}/${TENSOR} -o ${O} --outname ${N} \
            -t ${TOYS_PER_JOB} --toysDataMode ${MODE} --seed $((SEED_BASE + k)) \
            ${PM} ${EXTRA} > ${O}/log_${N}.log 2>&1
        python3 compact_toys.py "${O}/${N}.hdf5" 1 1 ${TOYS_PER_JOB} \
            -o "${SUM}" --delete-originals >> ${O}/log_${N}.log 2>&1
    ) &
done
wait
echo "### ${NAME} done"

python3 - "${O}" "${NAME}" "${JOBS}" <<'PY'
import sys
import numpy as np
O, name, jobs = sys.argv[1], sys.argv[2], int(sys.argv[3])
q, ndf = [], []
missing = []
for k in range(1, jobs + 1):
    try:
        d = np.load(f"{O}/toyGoF_{name}_batch{k}_summary.npz")
        q.append(d["q"]); ndf.append(d["ndf"])
    except FileNotFoundError:
        missing.append(k)
if missing:
    print(f"WARNING: missing batches {missing}")
if q:
    q, ndf = np.concatenate(q), np.concatenate(ndf)
    np.savez(f"{O}/toyGoF_{name}_summary.npz", q=q, ndf=ndf)
    print(f"combined summary: {len(q)} toys -> {O}/toyGoF_{name}_summary.npz")
PY
