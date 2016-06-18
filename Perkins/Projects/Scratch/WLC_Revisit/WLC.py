# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from FitUtil.WormLikeChain.Python.Code import WLC_Fit
from GeneralUtil.python import PlotUtilities as pPlotUtil
def Key(i,j):
    return str(i) + str(j)

if __name__ == "__main__":
    base = "./Out/"
    # only Lp/kbT matters, so just set kbT=1
    kbT = 4.1e-21
    N = 10
    LpArr = np.linspace(1e-9,60e-9,num=N)
    L0Arr = np.linspace(0.7,1.1,num=N)
    K0 = 1000e-12
    ext = np.linspace(0,1.3,num=500)
    ForceArr = dict()
    for i,Lp in enumerate(LpArr):
        for j,L0 in enumerate(L0Arr):
            Force = WLC_Fit.WlcExtensible(ext,kbT,Lp,L0,K0)
            ForceArr[Key(i,j)] = Force.copy()
    # POST: have all the forces, plot
    for i,Lp in enumerate(LpArr):
        for j,L0 in enumerate(L0Arr):
            MyForce = ForceArr[Key(i,j)]
            fig = pPlotUtil.figure()
            plt.plot(ext,MyForce)
            pPlotUtil.lazyLabel("Extension","Force","")
            Id = "{:d}-{:d}_Lp={:.3g}_L0={:.3g}.png".format(i,j,Lp,L0)
            pPlotUtil.savefig(fig,base + Id)
