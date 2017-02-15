# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import interpolate
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from GeneralUtil.python import PlotUtilities
from scipy.stats import norm
import Analysis

def plot_distribution(x,f,f_interp,f_cdf,thresh,num_bins=50):
    idx_high = np.where(f_cdf >= thresh)
    idx_low = np.where(f_cdf <= thresh)
    events = f_cdf[idx_low]
    f_minus_mu = f - f_interp
    mu,std = Analysis.spline_residual_mean_and_stdev(f,f_interp)
    limits = [min(f_minus_mu),max(f_minus_mu)]
    num_bins = 50
    linspace_f_diff = np.linspace(*limits,num=num_bins)
    pdf_diff = norm.pdf(x=linspace_f_diff,loc=mu,scale=std)
    plt.subplot(3,1,1)
    plt.plot(x,f_interp,color='b',linewidth=3)
    plt.plot(x,f,color='k',alpha=0.3)
    PlotUtilities.lazyLabel("Time (au)","Force (au)","")
    plt.subplot(3,1,2)
    plt.hist(f_minus_mu,bins=num_bins,normed=True)
    plt.plot(linspace_f_diff,pdf_diff,linewidth=3,color='r')
    PlotUtilities.lazyLabel("Force Difference (au)","Probability (au)","")
    plt.subplot(3,1,3)
    plt.semilogy(x[idx_high],f_cdf[idx_high],alpha=0.3,color='k',
                 label="Non-Event")
    plt.semilogy(x[idx_low],f_cdf[idx_low],linestyle='None',marker='.',
                 color='r',label="Event")
    PlotUtilities.lazyLabel("Time (au)","Probability","",frameon=True,
                            loc="lower right")

def plot_autocorrelation_log(x,*args):
    """
    plots the autocorrelation function and fit
    
    Args:
        x: corrlation abscissa
        *args: output of auto_correlation_tau
    Returns:
        nothing, plots the autocorrelation log 
    """
    tau,coeffs,auto = args
    tol = 1e-6
    auto_norm = (auto-min(auto))/(max(auto)-min(auto))
    log_norm = np.log(auto_norm + tol)
    plt.plot(x,log_norm)
    plt.plot(x,np.polyval(coeffs,x=x))
    plt.ylim([np.percentile(log_norm,0.5),max(log_norm)])

    
    

