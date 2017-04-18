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

delete_checkpoints_in(){
    dir="$1"
    rm -f "$dir"/*folds*.pkl 
    rm -f "$dir"/*Scores.pkl
    rm -f "$dir"/*param*.pkl
}
# Description:
# removes all checkpoints from the cache directories

delete_checkpoints_in cache_protein
delete_checkpoints_in cache
rm -f debug_no_event/*
rm -f debug_no_event_protein/*



