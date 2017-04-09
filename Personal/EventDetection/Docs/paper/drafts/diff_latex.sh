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
# runs latex git diff with the revision we care about, ignoring errors (
# mostly due to figure problems)
# Arguments:
#### Arg 1: Description
#rev=d61409367a84cd7389ecc1cd97d37bd0369520d1
rev=0f0d13f088c4bba06a3566f4a5efbcda3b6ef161
#  need to use --append-context2cmd to tract abstract, see
# researchgate.net/post/How_can_I_track_changes_in_LaTeX_especially_abstractenvironment
context_append_str="abstract,appendix,bibliography"
# essentially, ignore lstlisting:
# tex.stackexchange.com/questions/136786/latexdiff-and-lstlistinga
config_str='PICTUREENV=(?:picture|DIFnomarkup|lstlisting)[\w\d*@]*'
git-latexdiff -b \
	      --ignore-latex-errors \
	      --main prheenan_thesis_ms.tex $rev \
	      -- \
	      --config=$config_str \
	      --append-context2cmd=$context_append_str

