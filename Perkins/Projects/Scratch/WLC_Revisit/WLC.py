# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from FitUtil.WormLikeChain.Python.Code import WLC_ComplexValued_Fit 


from GeneralUtil.python import PlotUtilities as pPlotUtil


def run():
    ForceInitial = np.linspace(0,250,1000) * 1e-12
    kbT = 4.1e-21
    Lp = 10e-9
    L0 = 650e-9
    K0 = 1000e-12
    LpFactors = [1,1.2,2,5,10,20,50,100]
    Extensions = []
    Forces = []
    for Factor in LpFactors:
        Args = (kbT,Lp*Factor,L0,K0,ForceInitial)
        ExtensionTmp = WLC_ComplexValued_Fit.ExtensionPerForce(*Args)
        FinalForce = WLC_ComplexValued_Fit.InvertedWlcForce(ExtensionTmp,*Args)
        Forces.append(FinalForce)
        Extensions.append(ExtensionTmp)
    # now plot everything
    XPlotFunc = lambda x: x*1e9
    YPlotFunc = lambda y: y*1e12
    for i,(Ext,Force) in enumerate(zip(Extensions,Forces)):
        plt.plot(XPlotFunc(Ext),YPlotFunc(Force),
                 label="Lp~{:.2f}".format(LpFactors[i]))
    pPlotUtil.lazyLabel("Extension [nm]","Force [pN]","")
    plt.show()

if __name__ == "__main__":
    run()
