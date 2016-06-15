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
    Utility file: we write down our stocks and desired concentrations
    """
    stocks = np.array([96,166])
    # essentially, how many tubes did we combine?
    volumes = np.array([68,33*4])
    # For loading prep gel, dont want more than 20ug=20e3ng in 125uL, giving
    # 160 ng/uL
    DesiredConc = 50
    obj = DilutionUtil.PrintDilutions(stocks,volumes,DesiredConc)

if __name__ == "__main__":
    run()
