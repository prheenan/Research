# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import interpolate
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

class split_force_extension:
    """
    class representing a force-extension curve, split into approach, dwell,
    and retract
    """
    def __init__(self,approach,dwell,retract):
        self.approach = approach
        self.dwell = dwell
        self.retract = retract
    def zero_all(self,separation,zsnsr,force):
        self.approach.offset(separation,zsnsr,force)
        self.dwell.offset(separation,zsnsr,force)
        self.retract.offset(separation,zsnsr,force)
    def flip_forces(self):
        self.approach.LowResData.force *= -1
        self.dwell.LowResData.force *= -1
        self.retract.LowResData.force *= -1



def zero_by_approach(split_fec,n_smooth,flip_force=True):
    """
    zeros out (raw) data, using n_smooth points to do so 
    
    Args:
        split_fec: instead of split_force_extension
        n_smooth: number of points for smoothing
        flip_force: if true, multiplies the zeroed force by -1
    Returns:
        nothing, but modifies split_fec to be zerod appropriately. 
    """
    # PRE: assume the approach is <50% artifact and invols
    approach = split_fec.approach
    force_baseline = np.median(approach.Force)
    filtered_approach = FEC_Util.GetFilteredForce(approach,n_smooth)
    # find the last place we are above the median
    idx_where = np.where(filtered_approach.Force <= force_baseline)[0]
    idx_surface = idx_where[-1]
    # iterate once in order to get a better estimate of the baseline; we can
    # remove the effect of the invols entirely 
    force_baseline = np.median(approach.Force[:idx_surface])
    idx_surface =  np.where(filtered_approach.Force <= force_baseline)[0][-1]
    # get the separation at the baseline
    separation_baseline = filtered_approach.Separation[idx_surface]
    zsnsr_baseline = filtered_approach.Zsnsr[idx_surface]
    # zero everything 
    split_fec.zero_all(separation_baseline,zsnsr_baseline,force_baseline)
    if (flip_force):
        split_fec.flip_forces()
   
        
def split_FEC_by_meta(time_sep_force_obj):
    """
    given a time_sep_force object, splits it into approach, retract, and dwell
    by the meta information
        
    Args:
        time_sep_force_obj: whatever object to split, should have triggertime
        and dwelltime 
    Returns:
        scipy.interpolate.LSQUnivariateSpline object, interpolating f on x
    """
    start_of_dwell_time = time_sep_force_obj.TriggerTime
    end_of_dwell_time = start_of_dwell_time + \
                        time_sep_force_obj.SurfaceDwellTime
    get_idx_at_time = lambda t: np.argmin(np.abs(time_sep_force_obj.Time-t))
    start_of_dwell = get_idx_at_time(start_of_dwell_time)
    end_of_dwell = get_idx_at_time(end_of_dwell_time)
    # slice the object into approach, retract, dwell
    slice_func = lambda s: \
        FEC_Util.MakeTimeSepForceFromSlice(time_sep_force_obj,s)
    approach = slice_func(slice(0             ,start_of_dwell,1))
    dwell    = slice_func(slice(start_of_dwell,end_of_dwell  ,1))
    retract  = slice_func(slice(end_of_dwell  ,None          ,1))
    return split_force_extension(approach,dwell,retract)

def spline_interpolator(tau_x,x,f,deg=2):
    """
    returns a spline interpolator with knots uniformly spaced at tau_x over x
    
    Args:
        tau_x: the step size in whatever units of c
        x: the unit of 'time'
        f: the function we want the autocorrelation of
        deg: the degree of the spline interpolator to use. continuous to 
        deg-1 derivative
    Returns:
        scipy.interpolate.LSQUnivariateSpline object, interpolating f on x
    """
    # note: stop is *not* included in the iterval, so we add add an extra strep
    step_knots = tau_x/2
    knots = np.arange(start=min(x),stop=max(x)+step_knots,step=step_knots)
    # get the spline of the data
    spline_args = \
        dict(
            # degree is k, (k-1)th derivative is continuous
            k=deg,
            # specify the spline knots (t) uniformly in time at the 
            # autocorrelation time. dont want the endpoints
            t=knots[1:-1]
            )
    return interpolate.LSQUnivariateSpline(x=x,y=f,**spline_args)

def auto_correlation_tau(x,f_user,deg_autocorrelation=1):
    """
    get the atucorrelation time of f (ie: fit polynomial to log(autocorrelation)
    vs x, so the tau is more or less the exponential decay constant
    
    Args:
        x: the unit of 'time'
        f_user: the function we want the autocorrelation of 
        deg_autocorrelation: the degree of autocorrelation to use. defaults to  
        linear, to get the 1/e time of autocorrelation
    Returns:
        tuple of <autocorrelation tau, coefficients of log(auto) vs x fit,
                  auto correlation ('raw')>
    """
    f = f_user - np.mean(f_user)
    auto = np.correlate(f,f,mode='full')
    # only want the last half (should be identical?) 
    size = int(auto.size/2)
    auto = auto[size:]
    # normalize the auto correlation, add in a small bias to avoid 
    # taking the log of 0. data is normalized to 0->1, so it should be OK
    tol = 1e-9
    # auto norm goes from 0 to 1
    auto_norm = (auto - np.min(auto))/(np.max(auto)-np.min(auto)) 
    # statistical norm goes from -1 to 1
    statistical_norm = (auto_norm - 0.5) * 2
    log_norm = np.log(auto_norm + tol)
    median = np.median(log_norm)
    fit_idx_max = np.where(statistical_norm <= 0)[0]
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
    return tau,coeffs,auto
    
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

    
    

