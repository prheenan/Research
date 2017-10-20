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

dir="//perknas2.colorado.edu/group/4Patrick/CuratedData/DNA/prc2-dna-binding/"
dir+="prc2-binding/10mM-MgCl2/analysis/in-0x/"

N=20
cd $dir
for x in `find $dir \( -iname "*.txt" -o -iname "*.pkl" -o -iname "*.tiff" \)`; do
    file_name=`basename $x`
    len=${#file_name}
    echo $file_name    
    offset_from_end=`expr $len - $N`
    name_start=`echo $file_name | cut -c -$N`
    name_id=`echo $file_name | grep -oP '(Image\d+)'`
    name_protein=`echo $file_name | grep -oP '(DNA|Protein_DNA)?.{0,4}?(txt|pkl|tiff)$'`
    new_name="$name_start-$name_id-$name_protein"
    echo "$file_name to $new_name"
    echo "___"
    mv -f $file_name $new_name || echo "Same"
done
cd - 
# Returns:



