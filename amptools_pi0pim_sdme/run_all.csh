#!/bin/tcsh -f
# One command for the whole chain:  ./run_all.csh my_sdme.ini [extra amptools_sdme options]
# Set up your GlueX/ROOT environment first (so `root` and the halld_sim `fit` work).

set here = "$0:h"
if ( "$here" == "$0" ) set here = .
pushd "$here" > /dev/null
set here = "$cwd"
popd > /dev/null
if ( $#argv >= 1 ) then
    set ini = "$1"
    shift
else
    set ini = example_sdme.ini
endif

python3 -c "import uproot, matplotlib" >& /dev/null
if ( $status != 0 ) then
    echo "Installing uproot and matplotlib ..."
    python3 -m pip install --user uproot matplotlib
    if ( $status != 0 ) exit 1
endif

python3 "$here/amptools_sdme.py" all -c "$ini" $argv:q
