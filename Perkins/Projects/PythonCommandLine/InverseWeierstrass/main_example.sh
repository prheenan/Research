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

set -x
python main_iwt.py \
    -number_of_pairs 16 \
    -number_of_bins  80 \
    -f_one_half 10e-12 \
    -fraction_velocity_fit 0.5 \
    -flip_forces 0 \
    -file_input ./Examples/input.pxp \
    -file_output landscape.csv



