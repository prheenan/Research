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
    # get the difference between what we expect and 
    # what we get
    f_minus_mu = f_interp - f
    limits = [min(f_minus_mu),max(f_minus_mu)]
    bins = 50
    linspace_f_diff = np.linspace(*limits,num=bins*5)
    # symetrically choose percentiles for the fit
    start_q = 1
    qr_1,qr_2 = np.percentile(a=f_minus_mu,q=[start_q,100-start_q])
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
    # make a threshold in probability (this will likely be machine-learned) 
    thresh = 1e-4
    idx_high = np.where(force_cdf >= thresh)
    idx_low = np.where(force_cdf <= thresh)
    events = force_cdf[idx_low]
    # plot everything
    fig = PlotUtilities.figure(figsize=(8,16))
    plt.subplot(3,1,1)
    plt.plot(x,f_interp,color='b',linewidth=3)
    plt.plot(x,f,color='k',alpha=0.3)
    PlotUtilities.lazyLabel("Time (au)","Force (au)","")
    plt.subplot(3,1,2)
    plt.hist(f_minus_mu,bins=bins,normed=True)
    plt.plot(linspace_f_diff,pdf_diff,linewidth=3,color='r')
    PlotUtilities.lazyLabel("Force Difference (au)","Probability (au)","")
    plt.subplot(3,1,3)
    plt.semilogy(x[idx_high],force_cdf[idx_high],alpha=0.3,color='k',
                 label="Non-Event")
    plt.semilogy(x[idx_low],force_cdf[idx_low],linestyle='None',marker='.',
                 color='r',label="Event")
    PlotUtilities.lazyLabel("Time (au)","Probability","",frameon=True,
                            loc="lower right")
    PlotUtilities.savefig(fig,"out.png")



if __name__ == "__main__":
    run()
