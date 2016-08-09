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

def GetPrimersWithoutOverhangs(Plasmid,ProductLength,
                               SliceOther,OtherIsReverse,PrimerLength,Overhang,
                               Name):
    # overhang is length + abasic site
    OverhangLength = len(Overhang) + 1
    # need zero based to get thi
    StartOther = SliceOther.start
    EndOther = SliceOther.stop
    LengthOtherPrimer = abs(EndOther-StartOther)+1
    TotalPrimerLength = PrimerLength +LengthOtherPrimer +OverhangLength
    # get the remaining length
    DistanceBetweenPrimers = ProductLength-TotalPrimerLength
    # figure out where *this* primer should start
    if (OtherIsReverse):
        # get the start
        StartOfOther = min(StartOther,EndOther)
        Start = StartOfOther - DistanceBetweenPrimers - PrimerLength
    else:
        Start = StartOther+ LengthOtherPrimer + DistanceBetweenPrimers
    End = Start + PrimerLength
    PrimerSlice = slice(Start,End,1)
    # determine which slice goes where
    if (OtherIsReverse):
        ReverseSlice = SliceOther
        ForwardSlice = PrimerSlice
    else:
        ReverseSlice = PrimerSlice
        ForwardSlice = SliceOther
    # cool, now just print out the sequences
    PrimerForwardSeq = Plasmid[ForwardSlice]
    PrimerReverseSeq = KmerUtil.ReverseComplement(Plasmid[ReverseSlice])
    return PrimerForwardSeq,PrimerReverseSeq,ForwardSlice,ReverseSlice

def CreatePrimer(Plasmid,ProductLength,
                 SliceOther,OtherIsReverse,PrimerLength,Overhang,Name):
    PrimerForwardSeq,PrimerReverseSeq,ForwardSlice,ReverseSlice = \
        GetPrimersWithoutOverhangs(Plasmid,ProductLength,
                                   SliceOther,OtherIsReverse,PrimerLength,
                                   Overhang,Name)
    print("Forward Slice is: {:s}".format(ForwardSlice))
    print("Reverse Slice is: {:s}".format(ReverseSlice))
    # for without the overhang, dont add *anything*
    OverhangUtil.ConcatAndSave("",baseDir="./",Name=Name,
                               ForwardSequence=PrimerForwardSeq,
                               ReverseSequence=PrimerReverseSeq,
                               addSpacer=False,addDbcoAndBio=False)
    OverhangUtil.ConcatAndSave(Overhang,baseDir="./",Name=Name,
                               ForwardSequence=PrimerForwardSeq,
                               ReverseSequence=PrimerReverseSeq,
                               addSpacer=True,addDbcoAndBio=False)
    OverhangUtil.ConcatAndSave(Overhang,baseDir="./",Name=Name,
                               ForwardSequence=PrimerForwardSeq,
                               ReverseSequence=PrimerReverseSeq,
                               addSpacer=True,addDbcoAndBio=True)


def run():
    """
    Using the overhang we had previous, constructs a new primer

    the 5' (1670F) is generally (1) more expensive and (2) we get less of it,
    so I will use the 3' (Biotin-TEG) to make this primer, which will be 
    much closer than 3520R
    """
    Plasmid = EmbossUtil.ReadSimpleSeqFile(inputFile)
    Overhang = KmerUtil.ReverseComplement("TAGGACCACTCT")
    SliceReverse = slice(3490,3520,1)
    SliceForward = slice(1606,1606+27,1)
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
