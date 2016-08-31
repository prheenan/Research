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
    stocks = np.array([5,4])
    # what volume are the stocks, in uL
    volumes = np.array([7,175*0.25])
    # what the post-dilution concenration is, ng/uL
    DesiredConc = np.array([0.2,1])
    obj = DilutionUtil.PrintDilutions(stocks,volumes,DesiredConc,
                                      UnitConc=["mg/mL","mM"],
                                      StockName=["Strept","TCEP "])

if __name__ == "__main__":
    run()
