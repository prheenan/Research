# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python.IgorUtil import SavitskyFilter

from Research.Personal.EventDetection.Docs.ToyGraphs import SimulationUtil
from Research.Personal.EventDetection.Util import Analysis


def run():
    """
    This shows how to get the 1/e autocorrelation time from a WLC-looking thing
    """
    x,dx,stretch_kwargs,f = SimulationUtil.\
            normalized_force_extension(snr=1e3,points=300)
    tau,linear_auto_coeffs,auto_norm = Analysis.auto_correlation_tau(x,f)
    num_points = np.round(tau/dx)
    smoothed = SavitskyFilter(inData=f,nSmooth=num_points)
    pred = np.exp(np.polyval(linear_auto_coeffs,x=x))
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    plt.plot(x,f,color='k',alpha=0.2,label="raw")
    plt.plot(x,smoothed,color='r',label="smoothed to autocorrelation time")
    PlotUtilities.\
        lazyLabel("","Force (au)",
                  "1/e Autocorrelation time by fitting to initial decay")
    plt.subplot(2,1,2)
    label = "Small-time fit, 1/e time is {:.2f} (N={:d} points)".\
            format(tau,int(np.round(tau/dx)))
    plt.plot(x,pred,label=label,linewidth=5,color='b',linestyle="--")
    plt.plot(x,auto_norm,label="autocorrelation",color='g')
    plt.ylim([min(auto_norm),max(auto_norm)])
    PlotUtilities.lazyLabel("Time (au)","Autocorrelation (au)","")
    PlotUtilities.savefig(fig,"autocorrelation.png")

if __name__ == "__main__":
    run()
