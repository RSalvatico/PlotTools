#!/bin/bash
# HTCondor GoF toy production -- large-N version of run_gof_toys_lean.sh, for
# when local 16-way parallelism is too slow (e.g. the 50k-toy default here).
#
#   cmsenv && source setup.sh && source setup_rabbit.sh
#   cd GoFStudy_orig_noFTS_lowess
#   ./submit_toys_condor.sh CR_observed rabbit/ourCR_postfitmean.hdf5 observed \
#       "--freezeParameters tt-vcb --unblind xsec_.*"
#   ./submit_toys_condor.sh SR_observed rabbit/ourSR_postfitmean.hdf5 observed \
#       "--unblind tt-vcb xsec_.*"
#   condor_q                                  # watch the cluster drain
#   ./submit_toys_condor.sh --collect CR_observed
#   ./submit_toys_condor.sh --collect SR_observed
#   python3 ../analysis/plotToyGoF.py --npz rabbit/toyGoF_CR_observed_summary.npz \
#       --observed rabbit/CR_observed.hdf5 -o plots/toyGoF_CR_observed.png
#
# Job output (condor dir, collected summaries) lands under O (default
# "rabbit"); override with e.g. O=rabbit_5FS for a different tensor/fit
# namespace (such as the 5FS-only ttbar datacard) so it never collides with
# another study's toy batches of the same CR_observed/SR_observed names:
#   O=rabbit_5FS ./submit_toys_condor.sh CR_observed rabbit_5FS/ourCR_postfitmean.hdf5 \
#       observed "--freezeParameters tt-vcb --unblind xsec_.*"
#   O=rabbit_5FS ./submit_toys_condor.sh --collect CR_observed
#
# Do NOT pass --paramModel in <extra flags> -- PM below already provides it;
# doubling it crashes rabbit's fitter with "Duplicate parameter names"
# (hit and fixed once already building run_all.sh, see git log).
#
# Storage safety, one level past run_gof_toys_lean.sh: each job compacts its
# own toy batch (nllvalreduced/ndfsat only, ~KB) on the WORKER's own local
# scratch and only ever transfers that tiny .npz + a log back over AFS -- the
# ~1.1 MB/toy raw fit-result hdf5 (a few hundred MB per job at these batch
# sizes) never touches AFS at all. See run_gof_toys_lean.sh's header for the
# AFS-quota crash this design avoids.
#
# TOYS=50000 JOBS=250 ./submit_toys_condor.sh <name> <tensor> <mode> <extra>
# -> 250 condor jobs x 200 toys each, ~85 min/job at ~25s/toy (6h runtime cap,
# edit +MaxRuntime below if a slower fit needs more headroom).

set -e
BASE=$(cd "$(dirname "$0")" && pwd)
cd "${BASE}"
O=${O:-rabbit}
PM="--paramModel Mu --paramModel analysis.rabbit_models.FreeNorm ttbb,ttbj,tt2b,ttcc,ttcj,tt2c,ttLF"

# --------------------------------------------------------------- collect mode
if [ "$1" = "--collect" ]; then
    NAME=$2
    [ -z "${NAME}" ] && { echo "usage: $0 --collect <name>" >&2; exit 1; }
    python3 - "${O}" "${NAME}" <<'PY'
import glob, sys
import numpy as np

O, name = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(f"{O}/toyGoF_{name}_batch*_summary.npz"))
if not files:
    raise SystemExit(f"no batch summaries found under {O}/toyGoF_{name}_batch*_summary.npz "
                      f"-- cluster still running, or wrong --name?")
q, ndf = [], []
for f in files:
    d = np.load(f)
    q.append(d["q"]); ndf.append(d["ndf"])
q, ndf = np.concatenate(q), np.concatenate(ndf)
np.savez(f"{O}/toyGoF_{name}_summary.npz", q=q, ndf=ndf)
print(f"combined {len(files)} batches, {len(q)} toys -> {O}/toyGoF_{name}_summary.npz")
PY
    exit 0
fi

# ---------------------------------------------------------------- submit mode
NAME=$1 TENSOR=$2 MODE=$3 EXTRA=$4
if [ -z "${NAME}" ] || [ -z "${TENSOR}" ] || [ -z "${MODE}" ]; then
    echo "usage: $0 <name> <tensor.hdf5> <expected|observed> <extra rabbit_fit.py flags>" >&2
    echo "       $0 --collect <name>" >&2
    exit 1
fi

TOYS=${TOYS:-50000}
JOBS=${JOBS:-250}
TOYS_PER_JOB=$(( (TOYS + JOBS - 1) / JOBS ))
SEED_BASE=${SEED_BASE:-70000}

CMSSW_SRC=$(cd "${BASE}/../.." && pwd)          # .../CMSSW_15_0_10/src
PT=$(cd "${BASE}/.." && pwd)                    # .../CMSSW_15_0_10/src/PlotTools
TENSOR_ABS=$(cd "$(dirname "${TENSOR}")" && pwd)/$(basename "${TENSOR}")
[ -f "${TENSOR_ABS}" ] || { echo "tensor not found: ${TENSOR_ABS}" >&2; exit 1; }

CONDOR_DIR=${BASE}/${O}/condor_${NAME}
mkdir -p "${CONDOR_DIR}"

EXE=${CONDOR_DIR}/run_batch.sh
cat > "${EXE}" <<EOF
#!/bin/bash
set -e
SCRATCH=\$(pwd)
K=\$1
cd ${CMSSW_SRC} && eval \`scramv1 runtime -sh\`
cd ${PT} && source setup.sh > /dev/null 2>&1
# 'getenv = True' carries RABBIT_ACTIVE in from the submit shell, which makes
# setup_rabbit.sh early-return without re-activating -- and 'scramv1 runtime'
# just above has already pushed CMSSW's python 3.9 ahead of the venv in PATH.
# rabbit_fit.py then runs under that python, which sees cvmfs tensorflow but not
# the venv-only tensorflow_probability. Force a real re-activation.
unset RABBIT_ACTIVE
source setup_rabbit.sh > /dev/null 2>&1
cd "\${SCRATCH}"
N=toyGoF_${NAME}_batch\${K}
rabbit_fit.py ${TENSOR_ABS} -o . --outname \${N} \\
    -t ${TOYS_PER_JOB} --toysDataMode ${MODE} --seed \$((${SEED_BASE} + K)) \\
    ${PM} ${EXTRA} > \${N}.log 2>&1
python3 ${PT}/GoFStudy_orig_noFTS_lowess/compact_toys.py \${N}.hdf5 \${K} \${K} ${TOYS_PER_JOB} \\
    -o \${N}_summary.npz --delete-originals >> \${N}.log 2>&1
EOF
chmod +x "${EXE}"

SUB=${CONDOR_DIR}/submit.sub
cat > "${SUB}" <<EOF
batch_name              = toyGoF_${NAME}
executable              = ${EXE}
arguments               = \$(Process)
transfer_executable     = True
+AccountingGroup        = "group_u_CMST3.all"
getenv                  = True
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_output_files   = toyGoF_${NAME}_batch\$(Process)_summary.npz,toyGoF_${NAME}_batch\$(Process).log
initialdir              = ${BASE}/${O}
output                  = ${CONDOR_DIR}/job.out.\$(Cluster).\$(Process)
error                   = ${CONDOR_DIR}/job.err.\$(Cluster).\$(Process)
log                     = ${CONDOR_DIR}/job.log.\$(Cluster)
+MaxRuntime             = 21600
RequestCpus             = 1
queue ${JOBS}
EOF

echo "### ${NAME}: ${JOBS} jobs x ${TOYS_PER_JOB} toys (${TOYS} total, seeds ${SEED_BASE}-$((SEED_BASE + JOBS - 1)))"
echo "### submit file: ${SUB}"
condor_submit "${SUB}"
echo "### monitor:  condor_q -batch"
echo "### results land in ${O}/toyGoF_${NAME}_batch*_summary.npz as each job finishes"
echo "### once the cluster drains: $0 --collect ${NAME}"
