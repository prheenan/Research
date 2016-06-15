# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
baseDir = "../../"
sys.path.append(baseDir)

import networkx


import PrimerDesign.OverHangingPrimers.OverhangGeneration as \
    OverhangGeneration
import Util.EmbossUtil as EmbossUtil
import Util.KmerUtil as KmerUtil
import Util.IdtUtil as IdtUtil
import Util.CommonPrimerUtil as CommonPrimerUtil

from itertools import product

def run():
    """
    Gets the best possible 15-mer with a melting temerature of about 35, 
    low self-dimerization
    """
    inputFile = baseDir + "PlasmidData/Plasmids/mp13_plasmid_plasmid_seq.txt"
    # get the smalles kmer with a decent number of kmers
    desiredPrimerLen = 12
    desiredMeltTemp=35
    # get the sequences, with biotin and dbco
    BaseDir="../../"
    OverhangGeneration.\
        CreateOverhangsFor1607F(inputFile,BaseDir,
                                desiredPrimerLen=desiredPrimerLen,
                                desiredMeltTemp=desiredMeltTemp,
                                Name="Ovh2.0",
                                MakeSpacerFile=True,
                                MakeLabelledFile=True,
                                CheckAfterChoosing=True,
                                ChooseFunc = lambda x: x[0])

if __name__ == "__main__":
    run()
