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
# runs a simulated unit test following Szabo, 2010, and a 'real' data test

input_file=${1:-"./UnitTest/input.pxp"}
output_file=${2:-"./UnitTest/unit_test.csv"}
expected_file=${3:-"./UnitTest/expected_landscape.csv"}

# # Run the Szabo test
szabo_base="../../../../../../FitUtil/EnergyLandscapes/InverseWeierstrass/Python"
szabo_path="${szabo_base}/TestExamples/Testing/"
cd "$szabo_path"
ls *.py | tail -1 | xargs python  || { echo "Szabo test failed" ; exit; }
echo "===Szabo Test Passed==="
cd - 
# POST: Szabo test passed,  in original directory
# go one directory up (where the python file lives)
cd ..
# remove the output file, in case we already ran
rm -f "$output_file"
python main_iwt.py \
    -number_of_pairs 16\
    -flip_forces 0\
    -number_of_bins 150\
    -f_one_half 8e-12\
    -fraction_velocity_fit 0.1\
    -file_input "${input_file}"\
    -file_output "${output_file}" || { echo "===IWT run failed===" ; exit; }
echo "===Run Test Passed ==="    
# make sure the output is what we expect
python ./UnitTest/test_equality.py "${output_file}" "${expected_file}" || \
    { echo "===Value test failed!==="; exit; }
echo "===Value Test Passed==="
                
