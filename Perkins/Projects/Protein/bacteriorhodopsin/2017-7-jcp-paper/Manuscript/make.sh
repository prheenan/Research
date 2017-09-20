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

make_figure()
{
    # see (or just man pandoc): pandoc.org/MANUAL.html
    # --filter      : use citeproc to generate citations, figure numbes
    # --bibliography: path to the zotero bibliography
    # --csl         : the style sheet to use with citeproc
    # --template    : the template to use 
    # --reference-docx: for getting better formatting
    # --from        : the type of format to use
    # --verbose     : print debugging info
    # -s            : make the file standalone
    # -o            : output
    # --metadata    : sets a relevant variable
    # note: metadata can be set as follows: 
    # stackoverflow.com/questions/26431719/pandoc-citations-without-appending-the-references-bibliography
    input_file="2017-jcp-bacteriorhodopsin.md" 
    pandoc $input_file style.yaml\
	--filter=./walk_figures.py\
	--filter=./fulfill_figures.py\
        --filter=pandoc-citeproc\
	--bibliography=${5:-./jcp_bR_energy_bibliography.bib}\
	--from=markdown+yaml_metadata_block\
	--csl=${4:-cite_style.csl}\
	--reference-docx=${3:-./Template/prh_template.docx}\
        --metadata link-citations=true\
	--verbose \
	-s -o "$input_file".docx
}

make_figure

# Arguments:
#### Arg 1: Description

# Returns:



