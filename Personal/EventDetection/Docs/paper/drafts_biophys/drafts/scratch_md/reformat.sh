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

input="supplemental.md"
output="${input}_sed.md"
output_tmp="${input}_tmp.md"
# note sure why, but I have to use strick alternation to get this to work...
# replace things like \citePRH{<XXX>} with just @
sed -E 's/\\cite[a-zA-Z]*\{([^}]+)\}/@\1/g' "$input" > $output
# replace ValUnit{<x>}{<y>} with '<x> <y>', same with supply<x><y>
sed -E 's/\\ValUnit{([^}]+)}{([^}]+)}/\1 \2/g' "$output" > $output_tmp
sed -E 's/\\supply{([^}]+)}{([^}]+)}/\1 \2/g' "$output_tmp" > $output
# replace eqlab{xx}{Label} with markdown equiv
sed -E 's/\\eqlab{([^}]+)}{([^}]+)}/$\1$ {#eq:\2}/g' "$output" > $output_tmp
# replace \chapter<x> and \section<x> with #x and ## x
sed -E 's/\\sLabel{([^}]+)}/ {#sec:\1 /g' "$output_tmp" > $output
sed -E 's/\\chapter{([^}]+)}/#\1/g' "$output" > $output_tmp
sed -E 's/\\section{([^}]+)}/##\1/g' "$output_tmp" > $output

mv $output $output_tmp

# replace \singlemol{} with SMFS
sed -iE "s/singlemol{}/SMFS/g" $output_tmp
# replace \name{} with FEATHER
sed -iE "s/name{}/FEATHER/g" $output_tmp
# replace \tRef/sRef/eRef/fRef with the equiv
sed -iE 's/\\tRef{/{#tbl:/g' $output_tmp
sed -iE 's/\\sRef{/{#sec:/g' $output_tmp
sed -iE 's/\\fRef{/{#fig:/g' $output_tmp
sed -iE 's/\\eRef{/{#eq:/g' $output_tmp
# replace degree C with whatever...
sed -iE 's/\\degreeC{}/$^{\\circ}$C/g' $output_tmp

# replace \ValUnit{}{} with just the actual words
# remove all the backslashes
mv $output_tmp $output


cat $output
#gsed -iE "s/chapter\{([^}]+)\}/@\1/g" $output


