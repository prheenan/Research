# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

baseDir = "../../"
sys.path.append(baseDir)
inputFile = baseDir + "PlasmidData/Plasmids/mp13_plasmid_plasmid_seq.txt"
import Util.EmbossUtil as EmbossUtil
import Util.CommonPrimerUtil as CommonPrimerUtil
import Util.KmerUtil as KmerUtil
from PrimerDesign.OverHangingPrimers.OverhangGeneration \
    import CreatePrimer


def run():
    """
    Using the overhang we had previous, constructs a new primer

    the 5' (1670F) is generally (1) more expensive and (2) we get less of it,
    so I will use the 3' (Biotin-TEG) to make this primer, which will be 
    much closer than 3520R
    """
    Plasmid = EmbossUtil.ReadSimpleSeqFile(inputFile)
    Overhang = KmerUtil.ReverseComplement("TAGGACCACTCT")
    SliceReverse = CommonPrimerUtil.SliceReverse
    SliceForward = CommonPrimerUtil.SliceForward
    ProductLength = 201
    """
    The following shows 200 should work well, *including* the nick/abasic site
    (j factor is a mechanical property)

    Vafabakhsh R, Ha T. 
    Extreme Bendability of DNA Less than 100 Base Pairs Long 
    Revealed by Single-Molecule Cyclization. 
    Science 2012
    """
    CreatePrimer(Plasmid,ProductLength,
                 SliceOther=SliceForward,OtherIsReverse=False,PrimerLength=22,
                 Overhang=Overhang,Name="ovh2.7_200nt")
    CreatePrimer(Plasmid,ProductLength,
                 SliceOther=SliceReverse,OtherIsReverse=True,PrimerLength=22,
                 Overhang=Overhang,Name="ovh2.5_200nt")


if __name__ == "__main__":
    run()
