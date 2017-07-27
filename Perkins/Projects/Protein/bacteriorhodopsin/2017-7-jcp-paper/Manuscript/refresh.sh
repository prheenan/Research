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
# Re-calculates all the figures

# run each python file
figure_base_dir="../Figures/"
out_dir="./Figures/"
python_files=`find $figure_base_dir -name "*.py"`
for f in $python_files
do
    dir_name_f=`dirname $f`
    base_name_f=`basename $f`
    cd "$dir_name_f"
    python "$base_name_f"
    # change back
    cd -
done
# copy all the (newly made) figures 
find "$figure_base_dir" -name "*.png" -exec cp {} "$out_dir" \;
# re-make the docx file. 




