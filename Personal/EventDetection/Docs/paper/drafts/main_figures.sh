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
    png=${3}
    cd ${1}
    for i in *.py; do
	# try to run the pdf; ignore errors
	python2 "$i" || true
    done
    abs=$PWD
    for i in *.pdf; do
	cp "$i" ${2}
    done
    if [ $3 -eq 1 ] 
	then
	for i in *.pdf; do
	    out="${2}/_$i.png"
	    inkscape "${abs}/$i" -z --export-dpi=100 --export-area-drawing \
		--export-png="${out}" 
	done
    fi
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
pngs=0
base_dir_rel="../Figures/"
cartoon_dir="${base_dir_rel}FigureCartoon/"
timing_dir="${base_dir_rel}FigureTiming/"
prep_dir="${base_dir_rel}FigurePrep/"
rupture_dir="${base_dir_rel}FigureRupture/"
pres_dir="${base_dir_rel}Presentation/"
prev_dir="${base_dir_rel}FiguresPreviousWork/"
out_previous="../Finals_Presentation/"
# copy the metric table here and to the output, for the paper/pres to use
cp ${base_dir_rel}FigurePerformance_CS/*.tex .
cp metric_table.tex "$out_dir"
# copt the previous results/other relevant work paper 
copy_pdfs "$prev_dir" "$out_previous" $pngs
# copy the presenatation figures
copy_pdfs "${pres_dir}bhattacharya/"  $out_path $pngs
copy_pdfs "${pres_dir}noise_distribution/" $out_path $pngs
copy_pdfs "${pres_dir}domain-specific-path/" $out_path $pngs 
# copy the paper figures
copy_pdfs "${base_dir_rel}FigurePerformance_FullSet_FEATHER/" $out_path $pngs
copy_pdfs "${base_dir_rel}FigurePerformance_DistanceOnly/" $out_path $pngs
copy_pdfs "${base_dir_rel}FigurePerformance_per_algorithm/" $out_path $pngs
copy_pdfs "${prep_dir}" $out_path $pngs
copy_pdfs "${base_dir_rel}FigureAlgorithm/" $out_path $pngs
copy_pdfs "${base_dir_rel}FigureTuning/" $out_path $pngs
copy_pdfs "${timing_dir}" $out_path $pngs
copy_pdfs "${cartoon_dir}" $out_path $pngs
copy_pdfs "${rupture_dir}" $out_path $pngs
copy_pdfs "${base_dir_rel}FigurePerformance_CS/" $out_path $pngs




