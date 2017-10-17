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
files=`find ../ -name "*.py"`
for f in $files
    do
    file_name=`basename $f`
    file_dir=`dirname $f`
    cd $file_dir
    python $file_name || \
        { echo "Something went wrong with $file_name" ; exit 1; }
    cd - 
done

# Returns:



