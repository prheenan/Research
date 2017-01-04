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
    print("First, dilute stocks to get 'x' concentration")
    stock_conc = 10
    desired_conc = [1,0.1,0.01,0.001]
    volumes = [20]
    DilutionUtil.PrintSerialSteps(stock_conc,volumes,sorted(desired_conc)[::-1],
                                  ConcString="uM",
                                  BufferString="TE")
    print("Next, get gel lanes (add 2uL loading buffer to each lane)")
    # to get the highest concentration we want, dilute hairpin 1000x. 
    stock_conc = 1
    desired_conc = [0.1,0.01,0.001,0.0001]
    volumes = [10]
    DilutionUtil.PrintSerialSteps(stock_conc,volumes,sorted(desired_conc)[::-1],
                                  ConcString="x Hairpin",
                                  BufferString="x Comp")


if __name__ == "__main__":
    run()
