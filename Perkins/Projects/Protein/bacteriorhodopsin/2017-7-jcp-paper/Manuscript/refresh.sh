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

function re_read()
{
    python_files=`find $1 -name "*.py"`
    for f in $python_files
    do
        dir_name_f=`dirname $f`
        base_name_f=`basename $f`
        echo $dir_name_f
        cd "$dir_name_f"
        if [[ $f == *"Generate"* ]]; then
            echo "Ignoring ${f}"
        else
            python "$base_name_f"
        fi
       
        # change back
        cd -
    done
}

# Description:
# Re-calculates all the figures

# run each python file
force_re_run=0
figure_base_dir="../Figures/"
out_dir="./Figures/"
if [[ $force_re_run == 1 ]]; then
    re_read "$figure_base_dir"
fi
# copy all the (newly made) figures 
find "$figure_base_dir" -name "*.png" -exec cp {} "$out_dir" \;
# re-make the docx file. 




