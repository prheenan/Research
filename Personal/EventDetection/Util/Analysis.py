# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

def auto_correlation_tau(x,f,deg_autocorrelation=1):
    """
    XXX switch to getting in units of number of points.
    """
    auto = np.correlate(f,f,mode='full')
    # only want the last half (should be identical?) 
    size = int(auto.size/2)
    auto = auto[size:]
    # normalize the auto correlation, add in a small bias to avoid 
    # taking the log of 0. data is normalized to 0->1, so it should be OK
    tol = 1e-9
    auto_norm = (auto - np.min(auto))/(np.max(auto)-np.min(auto)) + tol
    log_norm = np.log(auto_norm)
    median = np.median(log_norm)
    fit_idx_max = np.where(log_norm <= median)[0]
    assert fit_idx_max.size > 0 , "autocorrelation doesnt dip under median"
    # get the first time we cross under the median
    fit_idx_max =  fit_idx_max[0]
    # git a high-order polynomial to the auto correlation spectrum, get the 1/e
    # time.
    coeffs = np.polyfit(x=x[:fit_idx_max],y=log_norm[:fit_idx_max],
                        deg=deg_autocorrelation)
    # get just the linear and offset
    linear_auto_coeffs = coeffs[-2:]
    # get tau (coefficient in the exponent, y=A*exp(B*t), so tau=1/B
    # take the absolute value, since tau is a decay, has a minus 
    tau = abs(1/linear_auto_coeffs[0])
    return tau,linear_auto_coeffs,auto_norm

