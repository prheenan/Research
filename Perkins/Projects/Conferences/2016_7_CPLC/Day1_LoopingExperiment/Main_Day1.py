# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


def run():
    """

    """
    Numer3 = [19,185,282,261,340,315,301,276,290,245]
    Denom3 = [642,473,538,428,450,424,381,365,371,303]
    off = 3
    Times3 = [0,off+0,off+4,off+8,off+17,off+23,off+31,off+35,off+44,off+53]

    Numer4 = [21,175,349,435,470,495,552,536,504,524,455]
    Denom4 = [922,811,906,869,809,728,757,696,663,674,639]
    offFour = off
    Times4 = [0,offFour+0,offFour+3,offFour+8,offFour+11,offFour+20,
              offFour+25,offFour+33,offFour+37,offFour+47,offFour+59]
    # get the population ratio
    Ratio3 = np.array(Numer3)/np.array(Denom3)
    Ratio4 = np.array(Numer4)/np.array(Denom4)
    fig = plt.figure()
    plt.plot(Times3,Ratio3,'ro-',label="Sample 3")
    plt.plot(Times4,Ratio4,'bs--',label="Sample 4")
    print("Times for Sample 3 in minutes: {:s}".format(Times3))
    print("Ratios for Sample 3: {:s}".format(Ratio3))
    print("Times for Sample 4 in minutes: {:s}".format(Times4))
    print("Ratios for Sample 4: {:s}".format(Ratio4))
    plt.legend(loc='upper left')
    plt.xlabel("Time (minutes)")
    plt.ylabel("Population in High fret [0,1]")
    plt.show()

if __name__ == "__main__":
    run()
