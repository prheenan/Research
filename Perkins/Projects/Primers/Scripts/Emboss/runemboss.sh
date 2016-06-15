#!/bin/bash
# courtesy of: http://redsymbol.net/articles/unofficial-bash-strict-mode/
# (helps with debugging)
# set -e: immediately exit if we find a non zero
# set -u: undefined references cause errors
# set -o: single error causes full pipeline failure.
set -euo pipefail
IFS=$'\n\t'

## This script runs emboss on a input file, getting primers within the specified
## melting temperature and length recquirements. It tries to find primers
## with low self-alignment scores, but it doesnt use IDTs method of local
## alignment

# Args:
# 1: inFile: fasta file we read in
# 2: outFile: where to put the resulting primers
# 3: primer size: (optional): how long the primer should be. Default: 12
# 4: minTm (optional): minimum melting temeprature. Default: 37
# 5: maxTm (optional): maximum melting temeprature. Default: 43
# 6: optTm (optional): optimal melting temperature. Default: 40



# datestring, used in many different places...
dateStr=`date +%Y-%m-%d:%H:%M:%S`
# see the emboss documentation:
# http://emboss.sourceforge.net/apps/cvs/emboss/apps/eprimer3.html
inFile=$1
outFile=$2
#   -minsize            integer    [18] Minimum acceptable length of a primer.
#                                  Must be greater than 0 and less than or
#                                  equal to MAX-SIZE. (Integer 1 or more)
#*  -maxsize            integer    [27] Maximum acceptable length (in bases) of
#                                  a primer. Currently this parameter cannot
#                                  be larger than 35. This limit is governed by
#                                  the maximum oligo size for which
#                                  Eprimer32's melting-temperature is valid.
#                                  (Integer up to 35)
# *  -optsize            integer    [20] Optimum length (in bases) of a primer
#                                  oligo. Eprimer32 will attempt to pick
#                                  primers close to this length. (Integer 0 or
#                                  more)
primerSize=${3:-12}
#*  -psizeopt           integer    [200] The optimum size for the PCR product.
#                                  0 indicates that there is no optimum product
#                                  size. (Integer 0 or more)
optPcrProductSize=0
#*  -dnaconc            float      [50.0] The nanomolar concentration of
#                                  annealing oligos in the PCR.	  
dnaConc=50 #nM
saltConc=50 #mM
# temperature ranges in celcius
#*  -mintm              float      [57.0] Minimum acceptable melting
#                                  temperature(Celsius) for a primer oligo.
#                                  (Any numeric value)
#*  -maxtm              float      [63.0] Maximum acceptable melting
#                                  temperature(Celsius) for a primer oligo.
#                                  (Any numeric value)
#*  -maxdifftm          float      [100.0] Maximum acceptable (unsigned)
#                                  difference between the melting temperatures
#                                 of the forward and reverse primers. (Any
#                                 numeric value)
minTm=${4:-37}
maxTm=${5:-43}
optTm=${6:-40}
#*  -task               menu       [1] Tell Eprimer32 what task to perform.
#                                  Legal values are 1: 'Pick PCR primers', 2:
#                                  'Pick forward primer only', 3: 'Pick reverse
#                                  primer only', 4: 'No primers needed'.
#                                  (Values: 1 (Pick PCR primers); 2 (Pick
#                                  forward primer only); 3 (Pick reverse primer
#                                  only); 4 (No primers needed))
task=2
diffTm=40 # between forward and reverse primer, usually dont care for ovh
maxRet=2000
#*  -prange             range      [100-300] The associated values specify the
#                                  lengths of the product that the user wants
#                                  the primers to create, and is a space
#                                  separated list of elements of the form
#                                  (x)-(y)
minLen=0
maxLen=7249 # used to bound the product size; doesn't need to be exact
# maximum number of repeats
#*  -maxpolyx           integer    [5] The maximum allowable length of a
#                                  mononucleotide repeat in a primer, for
#                                  example AAAAAA. (Integer 0 or more)
maxRepeat=3
#  -gcclamp            integer    [0] Require the specified number of
#                                  consecutive Gs and Cs at the 3' end of both
#                                  the forward and reverse primer.
#*  -ogcopt             float      [50.0] Internal oligo optimum GC percent.
#                                  (Any numeric value)
#*  -ogcmin             float      [20.0] Minimum allowable percentage of Gs
#                                  and Cs in an internal oligo. (Any numeric
#                                  value)
#*  -ogcmax             float      [80.0] Maximum allowable percentage of Gs
#                                  and Cs in any internal oligo generated by
#                                  Primer. (Any numeric value)
maxGC=75 # percentage GC content
minGC=0
maxSize=${primerSize}
#
#   -scorrection        menu       [1] Specifies the salt correction formula
#                                  for the melting temperature calculation.
#                                  (Values: 0 (Schildkraut & Lifson); 1
#                                  (SantaLucia); 2 (Owczarzy et al))
scorrection=1
#'''
#   -tmformula          menu       [1] Specifies details of melting temperature
#                                  calculation. (Values: 0 (Breslauer et al);
#                                  1 (SantaLucia))
#'''
formula=1
# lowest possible alignment score. used  for...
#
#   -selfend            float      [3.00] The maximum allowable 3'-anchored
#                                  global alignment score
#  -selfany            float      [8.00] The maximum allowable local alignment
#                                  score when testing a single primer for
#                                  (local) self-complementarity 
# note: omishbymax is from 'plasmidFile'
alignScore=0.0
gcclamp=1 #number of GCs at 3' end..
set -x
eprimer32 -auto -outfile=${outFile} \
	  -sequence ${inFile} -task ${task}\
	  -gcclamp ${gcclamp} \
	  -ogcmin ${minGC} -maxgc ${maxGC} -mingc ${minGC} \
	  -minsize ${primerSize} \
	  -maxsize ${primerSize} \
	  -optsize ${primerSize}  \
	  -psizeopt ${optPcrProductSize} -prange 12-$maxLen \
	  -mintm ${minTm} -maxtm ${maxTm} -opttm ${optTm} \
	  -selfend ${alignScore} -selfany ${alignScore} \
	  -oendself ${alignScore} -oanyself ${alignScore} \
	  -saltConc ${saltConc} \
	  -maxdifftm ${diffTm} -numreturn ${maxRet} \
	  -maxpolyx ${maxRepeat} \
	  -tmformula ${formula} -scorrection ${scorrection} \
	  -dnaconc ${dnaConc} 
