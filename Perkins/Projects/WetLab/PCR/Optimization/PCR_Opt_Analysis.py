# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from Util import DilutionUtil

styles = [dict(color='r',marker='*',linestyle='--'),
          dict(color='b',marker='o',linestyle='-'),
          dict(color='g',marker='v',linestyle='-.'),
          dict(color='k',marker='s',linestyle='-',linewidth=3.0),
          dict(color='m',marker='x',linestyle='--',linewidth=3.0),
          dict(color='y',marker='.',linestyle='-',linewidth=3.0)]

def PCR_Analyze_Objects(Objs,SaveName,PrintText=True,DesiredConc=50):
    """
    Wrapper to analyze the objects.

    Args:
        Objs: list of Gradient objects
        SaveName: Base name for saving plots
    """
    fig = plt.figure()
    nStyles = len(styles)
    # combine things with the same repeat and machine
    for i,o in enumerate(Objs):
        YieldPerTube = o.GetYieldPer100uLTube()
        # want ug, everything in ng
        YToUnits = lambda y: y/1000
        label = "{:d}x Repeats, {:s}".format(int(o.Repeats),o.GetMachineName())
        plt.plot(o.Temperatures,YToUnits(YieldPerTube),label=label,
                 **(styles[i % nStyles]))
        if (PrintText):
            mStr = ""
            DilutionUtil.PrintDilutions(o.Concentrations,o.Volumes,
                                        DesiredConc)
    plt.xlabel("Annealing Temperature ('C)")
    plt.ylabel("Yield per 100uL tube (micrograms)")
    plt.legend()
    plt.title("PCR Optimization for " + SaveName)
    fig.savefig(SaveName + ".png")
