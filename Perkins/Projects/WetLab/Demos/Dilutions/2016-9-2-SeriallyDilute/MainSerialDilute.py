# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys



sys.path.append("../../../")
from Util import DilutionUtil 

def run():
    """
    Serially dilute something
    """
    ConcString = "ng/uL"
    VolString = "uL"
    # stock concentration: start off at 100uM
    Stock = 760e3/346
    # Desired concentrations
    Desired = [Stock/10,Stock/100,Stock/1000]
    # desired volumes (for each)
    Volumes = [105,80,80]
    DilutionUtil.PrintSerialSteps(Stock,Volumes,sorted(Desired)[::-1],
                                  ConcString="ng/uL",BufferString="TE")
    

if __name__ == "__main__":
    run()
