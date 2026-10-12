#!/bin/tcsh -f
# One command for the whole chain:  ./run_all.csh my_run.ini [extra amptools_auto options]
# Set up your GlueX/ROOT environment first (so $AMPTOOLS and root-config work).
# Builds extract_fit_results_etapi (needed by convert_to_csv.py) the first time,
# checks the Python packages the plotter needs, then runs every step.

set here = "$0:h"
if ( "$here" == "$0" ) set here = .
pushd "$here" > /dev/null
set here = "$cwd"
popd > /dev/null
if ( $#argv >= 1 ) then
    set ini = "$1"
    shift
else
    set ini = example_run.ini
endif

if ( ! -x "$here/extract_fit_results_etapi" ) then
    if ( ! $?AMPTOOLS ) then
        echo "AMPTOOLS is not set; set up your GlueX environment first."
        exit 1
    endif
    echo "Building extract_fit_results_etapi ..."
    g++ -O2 "$here/extract_fit_results_etapi.cc" -I"$AMPTOOLS" -L"$AMPTOOLS/lib" -lAmpTools \
        `root-config --cflags --libs` -o "$here/extract_fit_results_etapi"
    if ( $status != 0 ) exit 1
endif

python3 -c "import uproot, awkward, matplotlib" >& /dev/null
if ( $status != 0 ) then
    echo "Installing uproot, awkward and matplotlib for the plotter ..."
    python3 -m pip install --user uproot awkward matplotlib
    if ( $status != 0 ) exit 1
endif

python3 "$here/amptools_auto.py" all -c "$ini" $argv:q
