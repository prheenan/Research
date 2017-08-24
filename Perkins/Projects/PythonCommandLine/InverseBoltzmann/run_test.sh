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

# Returns:
python2 main_inverse_boltzmann.py \
    -number_of_bins 20\
    -interpolation_factor 1\
    -smart_interpolation 1\
    -gaussian_stdev 1e-8\
    -file_input Data/Experiment.pxp\
    -file_output ./out.csv



