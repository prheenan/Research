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

from scipy import interpolate
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
    tau,linear_auto_coeffs,auto_norm = Analysis.auto_correlation_tau(x,f)
    num_points = np.round(tau/dx)
    smoothed = SavitskyFilter(inData=f,nSmooth=num_points)
    # note: stop is *not* included in the iterval, so we add add an extra strep
    step_knots = tau/2
    knots = np.arange(start=min(x),stop=max(x)+step_knots,step=step_knots)
    # get the spline of the data
    spline_args = \
        dict(
            # degree is k, (k-1)th derivative is continuous
            k=2,
            # specify the spline knots (t) uniformly in time at the 
            # autocorrelation time. dont want the endpoints
            t=knots[1:-1]
            )
    interpolated = interpolate.LSQUnivariateSpline(x=x,y=f,**spline_args)
    f_interp = interpolated(x)
    f_deriv = interpolated.derivative()(x)
    ratio = -f_deriv/f_interp
    # get the difference between what we expect and 
    # what we get
    f_minus_mu = f_interp - f
    limits = [min(f_minus_mu),max(f_minus_mu)]
    bins = 50
    linspace_f_diff = np.linspace(*limits,num=bins*5)
    qr_1,qr_2 = np.percentile(a=f_minus_mu,q=[5,95])
    idx_fit = np.where( (f_minus_mu >= qr_1) &
                        (f_minus_mu <= qr_2))
    # fit a normal distribution to it, to get the standard deviation (globally)
    mu,std = norm.fit(f_minus_mu)
    mu,std = norm.fit(f_minus_mu[idx_fit])
    pdf_diff = norm.pdf(x=linspace_f_diff,loc=mu,scale=std)
    # get the distribution of the actual data
    distribution_force = norm(loc=f_interp, scale=std)
    # get the cdf of the data
    force_cdf = distribution_force.cdf(f)
    force_cdf_complement = 1-force_cdf
    # plot everything
    fig = PlotUtilities.figure()
    plt.subplot(3,1,1)
    plt.plot(x,f_interp,color='b',linewidth=3)
    plt.plot(x,f,color='k',alpha=0.3)
    PlotUtilities.lazyLabel("Time (au)","Force (au)","")
    plt.subplot(3,1,2)
    plt.hist(f_minus_mu,bins=bins,normed=True)
    plt.plot(linspace_f_diff,pdf_diff,linewidth=3,color='r')
    PlotUtilities.lazyLabel("Force Difference (au)","Probability (au)","")
    plt.subplot(3,1,3)
    plt.semilogy(x,force_cdf)
    ylim = [np.percentile(a=force_cdf,q=0.5),1]
    plt.ylim(ylim)
    PlotUtilities.lazyLabel("Time (au)","Probability","")
    PlotUtilities.savefig(fig,"out.png")



if __name__ == "__main__":
    run()
