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

function copy_pdfs()
{
    cd ${1}
    for i in *.py; do
	# try to run the pdf; ignore errors
	python2 "$i" || true
    done
    cd -
    cp ${1}*.pdf ${2}
}

# Description:

# Arguments:
#### Arg 1: Description

# Returns:

# get the output directory absolute (needed for inkscape)
out_dir="../Figures/Finals/"
mkdir -p $out_dir
cd $out_dir
out_path=$PWD
# go back and get the input directory absolute
cd - 
base_dir_rel="../Figures/"
cartoon_dir="${base_dir_rel}FigureCartoon/"
timing_dir="${base_dir_rel}FigureTiming/"
prep_dir="${base_dir_rel}FigurePrep/"
rupture_dir="${base_dir_rel}FigureRupture/"
#copy_pdfs ${base_dir_rel}FigurePerformance_FullSet_Only_FEATHER/ $out_path 
copy_pdfs ${prep_dir} $out_path 
copy_pdfs "${base_dir_rel}FigurePerformance_CS/" $out_path 
copy_pdfs "${base_dir_rel}FigureAlgorithm/" $out_path 
copy_pdfs "${base_dir_rel}FigureTuning/" $out_path 
copy_pdfs "${timing_dir}" $out_path 
copy_pdfs "${cartoon_dir}" $out_path 
copy_pdfs "${rupture_dir}" $out_path 



