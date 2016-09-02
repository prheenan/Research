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
    Stock = 43
    # Desired concentrations
    Desired = [2,1,0.5,0.2]
    # desired volumes (for each)
    Volumes = [25] + [10 for d in Desired[1:]]
    Stocks,VolumeStock,VolumeDilute,FinalStocks = DilutionUtil.\
                                    SeriallyDilute(Stock,Desired,Volumes)
    DilutionUtil.PrintSerialDilution(Stocks,VolumeStock,VolumeDilute,
                                     FinalStocks,ConcString=ConcString,
                                     VolString=VolString)
    

if __name__ == "__main__":
    run()
