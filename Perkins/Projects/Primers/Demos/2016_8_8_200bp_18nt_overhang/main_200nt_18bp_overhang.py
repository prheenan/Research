# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../")
import PrimerDesign.OverHangingPrimers.OverhangGeneration as Overhang
import Util.GraphUtil as GraphUtil
import Util.KmerUtil as KmerUtil
import Util.MeltingTemperatureUtil as MeltUtil
import Util.AlignUtil as AlignUtil
import Util.CommonPrimerUtil as CommonPrimerUtil
import networkx as nx


def run():
    """

    """
    InFile = "../../PlasmidData/Plasmids/mp13_plasmid_plasmid_seq.txt"
    OverhangLen = 18
    KmerLen = 6
    BaseDir = "../../"
    Choose = lambda x : x[0]
    SliceForward = CommonPrimerUtil.SliceForward
    ProductLength = 200
    inf,best = Overhang.CreateOverhangsFor1607F(InFile,BaseDir,OverhangLen,40,
                                                "Tmp",ChooseFunc=Choose,
                                                fudge=3,MakePrimerFile=False)
    Plasmid = inf.Plasmid
    # determine the minimum expected number of bases we would have in common
    NumBases = 4
    # this is how many bases we would expect to have in common for
    # a random plasmid
    NumBaseBits = np.ceil(np.log(len(Plasmid))/np.log(NumBases))
    # get the first (best) alignment score
    AlignmentScore = AlignUtil.GetEbiAlignments(str(best),Plasmid)[0].score
    print("Overhang has a plasmid alignment score of {:.2f}/{:.2f}/{:.2f}".\
          format(AlignmentScore,5*OverhangLen,NumBaseBits*5))
    PrimerLenAmplify = 28
    Name = "200nt_{:d}o_fwd".format(OverhangLen)
    Overhang.CreatePrimer(Plasmid,ProductLength,SliceOther=SliceForward,
                          OtherIsReverse=False,PrimerLength=PrimerLenAmplify,
                          Overhang=best,Name=Name)
                          
if __name__ == "__main__":
    run()
