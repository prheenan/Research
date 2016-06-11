# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


sys.path.append("../../")
from Util import DilutionUtil
from PCR.Optimization.PCR_Opt_Analysis import PCR_Repeat_Analysis


def run():
    """
    Demonstrates PCR optimization, given yields and gradient information
    """
    # in Celcius
    TemperaturesCommon = [60.0,61.5,62.3,64.0]
    Temperatures = [ TemperaturesCommon,
                     TemperaturesCommon]
    # number of repeats
    Repeats = [30,35]
    # Element Yields[i][j] is a list of ng/uL for Temperatures[j]
    # from Repeats[i]
    Yields = [ \
               [50.5,44.0,37.1,21.5],
               [85.1,85.9,66.7,48.3]
           ]
    Name = "Ovh4.0_Spacer"
    PCR_Repeat_Analysis(Temperatures,Repeats,Yields,Name,PrintText=True)
    


if __name__ == "__main__":
    run()
