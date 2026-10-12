#!/bin/bash
# One command for the whole chain:  ./run_all.sh my_run.ini [extra amptools_auto options]
# Set up your GlueX/ROOT environment first (so $AMPTOOLS and root-config work).
# Builds extract_fit_results_etapi (needed by convert_to_csv.py) the first time,
# checks the Python packages the plotter needs, then runs every step.
set -e
here=$(cd "$(dirname "$0")" && pwd)
ini=${1:-example_run.ini}; shift || true

if [ ! -x "$here/extract_fit_results_etapi" ]; then
    echo "Building extract_fit_results_etapi ..."
    g++ -O2 "$here/extract_fit_results_etapi.cc" -I"$AMPTOOLS" -L"$AMPTOOLS/lib" -lAmpTools \
        $(root-config --cflags --libs) -o "$here/extract_fit_results_etapi"
fi

if ! python3 -c "import uproot, awkward, matplotlib" 2>/dev/null; then
    echo "Installing uproot, awkward and matplotlib for the plotter ..."
    python3 -m pip install --user uproot awkward matplotlib
fi

python3 "$here/amptools_auto.py" all -c "$ini" "$@"
