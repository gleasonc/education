#!/bin/bash
# One command for the whole chain:  ./run_all.sh my_sdme.ini [extra amptools_sdme options]
# Set up your GlueX/ROOT environment first (so `root` and the halld_sim `fit` work).
set -e
here=$(cd "$(dirname "$0")" && pwd)
ini=${1:-example_sdme.ini}; shift || true

if ! python3 -c "import uproot, matplotlib" 2>/dev/null; then
    echo "Installing uproot and matplotlib ..."
    python3 -m pip install --user uproot matplotlib
fi

python3 "$here/amptools_sdme.py" all -c "$ini" "$@"
