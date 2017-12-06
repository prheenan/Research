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

# Description: updates the document (assuming it dosnt already exist) and the figures (always)

input_file="./2017-jcp-bacteriorhodopsin.md.docx"
output_base_dir="/c/Users/pahe3165/Dropbox/Perkins Group AFM/BR Energy landscape JPC/Additional Revision Documents/"
output_figure_dir="${output_base_dir}Figures/"
output_data_dir="${output_base_dir}Data/"
# cop the figures over
cp ./Figures/*.tiff "$output_figure_dir"
cp ./Figures/*.tiff "${output_figure_dir}tiff_only/"
cp ./Figures/*.svg "$output_figure_dir"
cp ./Figures/*.jpeg "$output_figure_dir"
cp ./Figures/*.csv "$output_data_dir"
