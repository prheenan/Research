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
from Research.Personal.EventDetection.Util import Analysis,Plotting

from scipy.stats import norm

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    x,dx,stretch_kwargs,f = SimulationUtil.\
            normalized_force_extension(snr=1e3,points=10000)
    tau,_,_ = Analysis.auto_correlation_tau(x,f)
    num_points = np.round(tau/dx)
    smoothed = SavitskyFilter(inData=f,nSmooth=num_points)
    interpolated = Analysis.spline_interpolator(tau,x,f)
    f_interp = interpolated(x)
    f_deriv = interpolated.derivative()(x)
    ratio = -f_deriv/f_interp
    # get the residual mean and standard deviation, from the spline...
    mu,std = Analysis.spline_residual_mean_and_stdev(f,f_interp)
    force_cdf = Analysis.spline_gaussian_cdf(f,f_interp,std)
    force_cdf_complement = 1-force_cdf
    # make a threshold in probability (this will likely be machine-learned) 
    thresh = 1e-4
    # plot everything
    fig = PlotUtilities.figure(figsize=(8,16))
    Plotting.plot_distribution(x,f,f_interp,force_cdf,thresh,num_bins=500)
    PlotUtilities.savefig(fig,"out.png")



if __name__ == "__main__":
    run()
