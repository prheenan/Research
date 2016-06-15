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
    Simply gets the 1607F and 3520R primers to-scale
    """
    # start off with a blank primer
    mPrimer = ""
    fwdSpacer,revSpacer =OverhangGeneration.\
                          ConcatPrimerTo1607Fand3520R(mPrimer,base=baseDir)
    # make orders, including the spaces
    mOrders = [[fwdSpacer,"1607F_4_2016"],
               [revSpacer,"3520R_4_2016"]]
    # must have 100nm scale for IdSp as of 4/1/2016
    scale = IdtUtil.Scales._100NM
    spacerOrders = IdtUtil.SequencesAndNamesTuplesToOrder(mOrders,
                                                          Scale=scale)
    IdtUtil.PrintAndSave(spacerOrders,"./IdtOrder_Linear.txt")
                                                           
if __name__ == "__main__":
    run()
