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
    OverhangGeneration.CreateOverhangsFor1607F(inputFile,baseDir,
                                               desiredPrimerLen=12,
                                               desiredMeltTemp=44,
                                               Name="Ovh3.0_L12_T50",
                                               MakeLabelledFile=True,
                                               MakeSpacerFile=True)
    # for the larger, want to build on our previous OVH 2.0 primer
    primer2p0 = "AGAGTGGTCCTA"
    # want to add an extra nucleotides to whatever 13-mer works
    # (note that the function checks we dont screw up when we return, so
    # this is okay. This is based on wanting a 5' GC clamp
    AddTo5Prime = "GA"
    AddTo3Prime = ""
    ChooseFuncLarger = lambda x: \
        OverhangGeneration.ChooseFirstPrimerMatching(x,primer2p0,
                                                     AddTo5Prime=AddTo5Prime,
                                                     AddTo3Prime=AddTo3Prime)
    OverhangGeneration.CreateOverhangsFor1607F(inputFile,baseDir,
                                               # note: we add 2, so total
                                               # length is 15
                                               desiredPrimerLen=13,
                                               desiredMeltTemp=44,
                                               # temperature with our add-on
                                               # (kludge) is about 45
                                               Name="Ovh4.0_L15_T52",
                                               MakeLabelledFile=True,
                                               ChooseFunc=ChooseFuncLarger,
                                               # we disable checking, since
                                               # we likely screw things up
                                               # with our addition
                                               CheckAfterChoosing=False,
                                               MakeSpacerFile=True)
    
if __name__ == "__main__":
    run()
