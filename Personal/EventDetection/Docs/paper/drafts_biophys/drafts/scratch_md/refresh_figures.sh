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

source ~/src_prh/GeneralUtil/bash/.profile
# Returns:
rm -f *feather*.docx
# copy all .png or .svg files in the protein directory
# (-o is OR clause)
find ../../../Figures/BiophysicalJournal/Protein/ \
    \( -name "*.png" -o -name "*.svg" \) -exec cp {} figures \;
# re-make the word document
p_pandoc scratch_prheenan.md supplemental_prheenan.md
date_str=`date +%Y_%m_%d-%H.%M`
out_name="feather_${date_str}.docx"
cp scratch_prheenan.md.docx $out_name
base_output="/Users/patrickheenan/Dropbox/Perkins Group AFM/FEATHER"
cp figures "$base_output/Figures/"
cp $out_name "$base_output/$out_name"
