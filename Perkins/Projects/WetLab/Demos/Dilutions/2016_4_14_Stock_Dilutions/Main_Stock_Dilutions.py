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
    stocks = [646]
    # what volume are the stocks, in <vol>
    volumes = [20]
    # what the post-dilution concenration is, <mass>/<vol>
    common_conc = 50
    DesiredConc = np.array([common_conc for _ in stocks])
    obj = DilutionUtil.PrintDilutions(stocks,volumes,DesiredConc,
                                      UnitConc=["nM" for _ in stocks])


if __name__ == "__main__":
    run()
