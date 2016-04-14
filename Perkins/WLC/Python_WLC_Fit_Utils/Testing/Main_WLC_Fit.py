# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

baseDir = "../"
sys.path.append(baseDir)
sys.path.append("../../../../../")
import WLC_Fit

def run():
    """
    Runs some unit testing on the WLC fitting
    """
    x = np.linspace(0,600e-9,15)
    xInterp = np.linspace(0,600e-9,100)
    y = np.array([0,0,0,0,0,0.5,1,1.5,3,4,9,30,50,100,250]) * 1e-12
    yInterp = np.interp(xInterp,x,y)
    mFit = WLC_Fit.ExtensibleWlcFit(xInterp,yInterp,VaryLp=False)
    mFitNon = WLC_Fit.NonExtensibleWlcFit(xInterp,yInterp,VaryLp=True)
    plt.plot(xInterp,yInterp,label="Extensible")
    plt.plot(xInterp,mFit,'r-',linewidth=3.0)
    plt.plot(xInterp,mFitNon,'b--',label="Non Extensible")
    plt.ylim([min(yInterp),max(yInterp)])
    plt.legend()
    plt.show()

if __name__ == "__main__":
    run()
