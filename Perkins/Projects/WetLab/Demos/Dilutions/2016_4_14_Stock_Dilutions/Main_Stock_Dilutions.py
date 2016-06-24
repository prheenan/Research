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
    stocks = np.array([13.2,54])
    # what volume are the stocks
    volumes = np.array([5,2])
    # what the post-dilution amount is 
    DesiredConc = 0.2
    obj = DilutionUtil.PrintDilutions(stocks,volumes,DesiredConc)

if __name__ == "__main__":
    run()
