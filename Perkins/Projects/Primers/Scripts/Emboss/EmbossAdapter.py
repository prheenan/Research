# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from subprocess import call

# XXX should fix this...
embossLoc= "../../Scripts/Emboss/runemboss.sh"

def RunEmbossAndCreatePrimerFile(inFile,outFile,length=12,minTm=33,maxTm=35,
                                 optTm=40):
    """
    Runs emboss on the specified input file with the given parameters, 
    saving out the primers to the outfile. All other parmaters are as defaulted
    in emboss util (e.g. gc clamp, etc)

    Args:
        inFile: where the input fasta file is
        outFile: where we should save the output primers
        length: length of the primer desired
        minTm: minimum melting temperature of the primer
        maxTm: maximum melting temperature of the primer
        optTm: optimal melting temperature of the primer 
    Returns:
        This is a description of what is returned.
    """
    args = [inFile,outFile,length,minTm,maxTm,optTm]
    strArgs = [str(a) for a in args]
    cmd =['bash',embossLoc] + strArgs
    call(cmd)
    # POST: file has been created, presummably. No need to return anything...
