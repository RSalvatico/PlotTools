#!/bin/bash
# Build rabbit_env/ from scratch. One-time, or after deleting a stale/broken env.
#
#   cmsenv
#   git submodule update --init rabbit   # if rabbit/ is empty
#   ./setup_rabbit_env.sh
#   source setup_rabbit.sh               # activate for day-to-day use
#
# CMSSW ships SciPy 1.10, whose minimizer callback signature silently disables
# the rabbit fit (every fit just returns its start values) -- pin scipy<1.14
# instead. After building, verify the env with a real fit whose answer differs
# from its start point (see COMBINE_FIXES.md).
#
# Never invoke the bare `pip` / `/usr/bin/pip` here: its shebang runs the system
# python, which under cmsenv's LD_LIBRARY_PATH loads CMSSW's mismatched
# libpython3.9 and dies with "No module named '_posixsubprocess'". CMSSW's
# `python3 -m venv` also builds the venv without pip. So we drive everything
# through "${VENV}/bin/python -m pip", which picks up CMSSW's pip 24.0 via
# --system-site-packages and still installs into the venv's own site-packages.

set -e
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="${BASE}/rabbit_env"
SRC="${BASE}/rabbit"
VPY="${VENV}/bin/python"

if [ ! -f "${SRC}/pyproject.toml" ]; then
    echo "setup_rabbit_env.sh: ${SRC} looks empty -- run 'git submodule update --init rabbit' first." >&2
    exit 1
fi

if [ -z "${CMSSW_BASE}" ]; then
    echo "setup_rabbit_env.sh: run 'cmsenv' first." >&2
    exit 1
fi

rm -rf "${VENV}"
python3 -m venv --system-site-packages "${VENV}"

# CMSSW's venv ships no pip; it must resolve via --system-site-packages.
if ! "${VPY}" -m pip --version; then
    echo "setup_rabbit_env.sh: '${VPY} -m pip' is not available -- is cmsenv active?" >&2
    exit 1
fi

"${VPY}" -m pip install -e "${SRC}[all]"
"${VPY}" -m pip install --ignore-installed --no-deps "scipy<1.14"

echo "rabbit_env built at ${VENV}. 'source setup_rabbit.sh' to activate, then run a real fit to verify."
