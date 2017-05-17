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
    # stock concentration
    Stock = 311
    # Desired concentrations
    Desired = [30,1,0.3,0.1]
    # desired volumes (for each)
    Volumes = [18,50,20,20,20]
    DilutionUtil.PrintSerialSteps(Stock,Volumes,sorted(Desired)[::-1],
                                  ConcString="ng/uL",BufferString="TE")
    

if __name__ == "__main__":
    run()
