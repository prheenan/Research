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
    Desired = [50,25,10]
    # desired volumes (for each)
    Volumes = [20 for d in Desired]
    DilutionUtil.PrintSerialSteps(Stock,Volumes,Desired,
                                  ConcString=ConcString,BufferString="TE",
                                  dilution_concentration=0)
    

if __name__ == "__main__":
    run()
