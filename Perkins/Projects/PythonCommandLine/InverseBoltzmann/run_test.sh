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

run()
{
    python2 main_inverse_boltzmann.py \
	-number_of_bins $1\
        -interpolation_factor 1\
        -output_interpolated $2\
        -smart_interpolation $3\
        -gaussian_stdev 1e-8\
        -file_input Data/Experiment.pxp\
        -file_output ./out.csv
}

# run several times, specifying the bins, smart interpolation, and 
run 20 1 1
run 20 1 0
run 20 0 1
run 20 0 0





