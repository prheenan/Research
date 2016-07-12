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
import Util.KmerUtil as KmerUtil
import PrimerDesign.OverHangingPrimers.OverhangGeneration as OverhangUtil


def run():
    """
    Using the overhang we had previous, constructs a new primer

    the 5' (1670F) is generally (1) more expensive and (2) we get less of it,
    so I will use the 3' (Biotin-TEG) to make this primer, which will be 
    much closer than 3520R
    """
    Plasmid = EmbossUtil.ReadSimpleSeqFile(inputFile)
    """
    The following shows 200 should work well, *including* the nick/abasic site
    (j factor is a mechanical property)

    Vafabakhsh R, Ha T. 
    Extreme Bendability of DNA Less than 100 Base Pairs Long 
    Revealed by Single-Molecule Cyclization. 
    Science 2012
    """
    DesiredLengthTotal = 200
    # overhang is 12nt + abasic site
    OverhangLength = 12 +1 
    # We want the new primer to have the same melting temperature as the
    # old one, so give it the same length. 5Primer starts (!) at 1607
    PrimerReverseLengthWithoutOverhang = 30
    # need zero based to get thi
    PrimerReverseStart = 3520
    PrimerReverseLengthWithOverhang = PrimerReverseLengthWithoutOverhang + \
                                      OverhangLength
    PrimerReverseSlice = slice(PrimerReverseStart-\
                               PrimerReverseLengthWithoutOverhang,
                               PrimerReverseStart,1)
    PrimerReverseSeq = KmerUtil.ReverseComplement(Plasmid[PrimerReverseSlice])
    # figure out where the 3' Primer should start
    PrimerForwardLengthWithoutOverhang = 23
    PrimerForwardLengthWithOverhang = PrimerForwardLengthWithoutOverhang + \
                                      OverhangLength
    LengthOfPrimersWithOverhangsAndAbasic = PrimerForwardLengthWithOverhang + \
                                            PrimerReverseLengthWithOverhang
    DistanceBetweenPrimers = DesiredLengthTotal - \
                             LengthOfPrimersWithOverhangsAndAbasic
    print(DistanceBetweenPrimers)
    PrimerForwardStart = PrimerReverseStart-DistanceBetweenPrimers-\
                         PrimerReverseLengthWithoutOverhang
    # okay, get the actual sequence we want
    PrimerSeq = slice(PrimerForwardStart-PrimerForwardLengthWithoutOverhang,
                      PrimerForwardStart,1)
    # Get the primer we want...
    SequencePrimerAnnealsTo = Plasmid[PrimerSeq]
    PrimerForward = SequencePrimerAnnealsTo
    # copied directly from IDT order on 2016/7/11, order from 2016-5-24
    Overhang = KmerUtil.ReverseComplement("TAGGACCACTCT")
    OverhangUtil.ConcatAndSave(Overhang,baseDir="./",Name="200nt_ovh2.5",
                               ForwardSequence=PrimerForward,
                               ReverseSequence=PrimerReverseSeq,
                               addSpacer=True,addDbcoAndBio=False)

if __name__ == "__main__":
    run()
