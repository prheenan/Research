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
    Stock = 20
    # Desired concentrations
    FirstDilutionConcentration = 10
    Desired = [FirstDilutionConcentration,2,0.4]
    # Maximum Volume possible to make
    StockVolumeTotal = 30
    MaxVolume = Stock * StockVolumeTotal / FirstDilutionConcentration
    # desired volumes (for each)
    Volumes = [20] + [20 for d in Desired[1:]]
    DilutionUtil.PrintSerialSteps(Stock,Volumes,Desired)
    

if __name__ == "__main__":
    run()
