# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../")
from Util import DilutionUtil 


def run():
    """
    fixed volume of hairpin, vary the relative concentration of DNA 
    """
    # to get the initial concentrations we want, need to dilute from 
    # 10uM to something like 1nM
    stock_conc = 100
    desired_conc = [10,1,1e-1,1e-2]
    volumes = [60]
    print("First, dilute hairpin stock to get 'x' concentration ({:.2f}uM)".\
          format(min(desired_conc)))
    DilutionUtil.PrintSerialSteps(stock_conc,volumes,sorted(desired_conc)[::-1],
                                  ConcString="uM",
                                  BufferString="TE")
    print("Next, get gel lanes (add 2uL loading buffer to each lane)")
    # to get all the concentrations we want, need 
    # (1) same volume of hairpin
    # (2) varying volumes of buffer and TE as appropriate...
    # note that this reduces the DNA concentration to about 90nM
    stock_conc_units_of_min = 1e4
    desired_conc = [1e3,3e2,1e2,3e1,1e1,1]
    volumes = [10]
    DilutionUtil.PrintSerialSteps(stock_conc_units_of_min,
                                  volumes,sorted(desired_conc)[::-1],
                                  ConcString="x Molar Comp",
                                  BufferString="1x W30R504T")


if __name__ == "__main__":
    run()
