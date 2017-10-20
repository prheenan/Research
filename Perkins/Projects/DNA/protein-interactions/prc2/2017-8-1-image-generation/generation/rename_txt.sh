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
dir+="prc2-binding/10mM-MgCl2/analysis/in_10mM"

dir="/c/Users/pahe3165/Desktop/tmp/"
N=20
cd $dir
ls
for x in `find $dir -name "*.txt"`; do
    file_name=`basename $x`
    len=${#file_name}
    echo $file_name    
    offset_from_end=`expr $len - $N`
    name_start=`echo $file_name | cut -c -$N`
    name_id=`echo $file_name | grep -oP 'Image\d+'`
    name_protein=`echo $file_name | grep -oP '(DNA|Protein_DNA)?\d*.{1,3}(txt|pkl)'`
    new_name="$name_start-$name_id$name_protein"
    echo $new_name
    echo "___"
    #mv $file_name $new_name
done
cd - 
# Returns:



