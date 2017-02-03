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
    stocks = [38]
    # what volume are the stocks, in <vol>
    volumes = [20]
    # what the post-dilution concenration is, <mass>/<vol>
    DesiredConc = np.array([38/5,38/15])
    obj = DilutionUtil.PrintDilutions(stocks,volumes,DesiredConc)


if __name__ == "__main__":
    run()
