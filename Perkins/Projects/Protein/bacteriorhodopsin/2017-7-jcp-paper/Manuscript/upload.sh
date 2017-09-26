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
output_base_dir="/Users/patrickheenan/Dropbox/Perkins Group AFM/BR Energy landscape JPC/"
output_figure_dir="${output_base_dir}Figures/"
output_dir="${output_base_dir}/drafts/"
output_file="Heenan_et_al_JCP_2017_bR_landscape_v1_prh.docx"
output_path="${output_dir}${output_file}"
# cop the figures over
cp ./Figures/*.png $output_figure_dir
if [ -f "$output_path" ]; then
    echo "File $output_path already exists. not-over-writing"
    exit
fi
cp "$input_file" "$output_path"
# Returns:



