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

function make_and_copy_figure()
{
    out_path=$1
    in_path_tmp=$2
    cd $in_path_tmp
    in_path=$PWD
    ext=".svg"
    # make all the PDFs!
    for i in *.svg; do
	in_tmp="$in_path/$i"
	# latex is unhappy with things like .svg.pdf, so just use .pdf
	in_without_extension=${i%.*}
	out_tmp="$out_path/$in_without_extension"
	# export_latex needed to avoid crappy rendering:
	# tex.stackexchange.com/questions/2099/how-to-include-svg-diagrams-in-latex
	inkscape "$in_tmp" --export-ps="${out_tmp}.ps"
	pstopdf "${out_tmp}.ps" -o "${out_tmp}.pdf"
	rm "${out_tmp}.ps"
    done
    cd -
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
make_and_copy_figure $out_path "${base_dir_rel}/FigureAlgorithm/"
make_and_copy_figure $out_path $timing_dir
make_and_copy_figure $out_path $cartoon_dir
make_and_copy_figure $out_path $prep_dir
make_and_copy_figure $out_path $rupture_dir
make_and_copy_figure $out_path "${base_dir_rel}FigurePerformance_CS/"



