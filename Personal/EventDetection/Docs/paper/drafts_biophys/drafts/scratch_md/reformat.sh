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

# Returns:

output=heenan_2017_ms_sed.tex 
# replace things like \citePRH{<XXX>} with just @
gsed -E "s/cite[a-zA-Z]*\{([^}]+)\}/@\1/g" heenan_2017_ms.tex > $output

# replace \singlemol{} with SMFS
gsed -iE "s/singlemol{}/SMFS/g" $output 
# replace \name{} with FEATHER
gsed -iE "s/name{}/FEATHER/g" $output 
# remove all the backslashes
gsed -i 's/\\//g' $output 
# replace \ValUnit{}{} with just the actual words
#gsed -ir "s/ValUnit\{([^}]+)\}\{([^}]+)\}/\1 \2/g" $output 


cat $output
#gsed -iE "s/chapter\{([^}]+)\}/@\1/g" $output


