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
output_tmp=heenan_2017_ms_sed_tmp.tex 
# replace things like \citePRH{<XXX>} with just @
sed -E "s/cite[a-zA-Z]*\{([^}]+)\}/@\1/g" heenan_2017_ms.tex > $output
sed -E "s/ValUnit{([^}]+)}{([^}]+)}/\1\2/g" $output > $output_tmp

# replace \singlemol{} with SMFS
sed -iE "s/singlemol{}/SMFS/g" $output_tmp
# replace \name{} with FEATHER
sed -iE "s/name{}/FEATHER/g" $output_tmp
# replace \ValUnit{}{} with just the actual words
# remove all the backslashes
sed -i 's/\/ /g' $output_tmp


cat $output
#gsed -iE "s/chapter\{([^}]+)\}/@\1/g" $output


