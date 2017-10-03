#!/bin/bash
# courtesy of: http://redsymbol.net/articles/unofficial-bash-strict-mode/
# (helps with debugging)
# set -e: immediately exit if we find a non zero
# set -u: undefined references cause errors
# set -o: single error causes full pipeline failure.
set -euo pipefail
IFS=$'\n\t'
# datestring, used in many different places...
dateStr=`date +%Y-%m-%d:%H:%M:%S`

# Description:

# Arguments:
#### Arg 1: Description
python_file="../../../PythonCommandLine/InverseWeierstrass/main_iwt.py"
python $python_file  \
    -number_of_pairs 10 \
    -number_of_bins 75 \
    -f_one_half 0 \
    -fraction_velocity_fit 0.1 \
    -flip_forces 1 \
    -file_input /Users/patrickheenan/src_prh/Research/Perkins/Projects/Trials/IWT/2016-7-26-devin/Hold.pxp \
    -file_output ./landscape.csv




