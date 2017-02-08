# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Personal.EventDetection.Docs.ToyGraphs import SimulationUtil


def run():
    """
    This shows how to get the 1/e autocorrelation time from a WLC-looking thing
    """
    x,dx,stretch_kwargs,f = SimulationUtil.normalized_force_extension(snr=1e3)
    auto = np.correlate(f,f,mode='full')
    # only want the last half (should be identical?) 
    size = int(auto.size/2)
    auto = auto[size:]
    # normalize the auto correlation, add in a small bias to avoid 
    # exponential problems
    tol = 1e-6
    auto_norm = (auto - np.min(auto))/(np.max(auto)-np.min(auto)) + tol
    log_norm = np.log(auto_norm)
    median = np.median(log_norm)
    fit_idx_max = np.where(log_norm <= median)[0]
    assert fit_idx_max.size > 0 , "autocorrelation doesnt dip under median?"
    # get the first time we cross under the median
    fit_idx_max =  fit_idx_max[0]
    # git a high-order polynomial to the auto correlation spectrum, get the 1/e
    # time.
    coeffs = np.polyfit(x=x[:fit_idx_max],y=log_norm[:fit_idx_max],deg=2)
    # get just the linear and offset
    linear_auto_coeffs = coeffs[-2:]
    pred = np.exp(np.polyval(linear_auto_coeffs,x=x))
    # get tau (coefficient in the exponent, y=A*exp(B*t), so tau=1/B
    # take the absolute value, since tau is a decay, has a minus 
    tau = abs(1/linear_auto_coeffs[0])
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    plt.plot(x,f)
    PlotUtilities.\
        lazyLabel("","Force (au)",
                  "1/e Autocorrelation time by fitting to initial decay")
    plt.subplot(2,1,2)
    label = "Small-time fit, 1/e time is {:.2f} (N={:d} points)".\
            format(tau,int(np.round(tau/dx)))
    plt.plot(x,pred,label=label)
    plt.plot(x,auto_norm,label="autocorrelation")
    plt.plot(x[fit_idx_max],pred[fit_idx_max],'ro')
    plt.ylim([min(auto_norm),max(auto_norm)])
    PlotUtilities.lazyLabel("Time (au)","Autocorrelation (au)",
                            "")
    PlotUtilities.savefig(fig,"autocorrelation.png")

if __name__ == "__main__":
    run()
