# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import interpolate
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from scipy.stats import norm

class split_force_extension:
    """
    class representing a force-extension curve, split into approach, dwell,
    and retract
    """
    def __init__(self,approach,dwell,retract,tau_num_points=None):
        self.approach = approach
        self.dwell = dwell
        self.retract = retract
        self.set_tau_num_points(tau_num_points)
    def retract_spline_interpolator(self,**kwargs):
        """
        returns an interpolator for force based on the stored time constant tau
        for the retract force versus time curbe

        Args:
            kwargs: passed to spline_interpolator
        """
        x,f = self.retract.Time,self.retract.Force
        return spline_interpolator(self.tau,x,f,**kwargs)
    def retract_separation_interpolator(self,**kwargs):
        """
        returns an interpolator for separation based on the stored time
        constant tau for the retract force versus time curbe

        Args:
            kwargs: passed to spline_interpolator
        """    
        x,f = self.retract.Time,self.retract.Separation
        return spline_interpolator(self.tau,x,f,**kwargs)
    def set_tau_num_points(self,tau_num_points):
        """
        sets the autocorrelation time associated with this curve
        
        Args:
            tau_num_points: integer number of points
        Returns:
            Nothing
        """
        self.tau_num_points = tau_num_points
        if (tau_num_points is not None):
            self.tau = np.median(np.diff(self.approach.Time))*tau_num_points
        else:
            self.tau = None
    def zero_all(self,separation,zsnsr,force):
        """ 
        zeros the distance and force of the approach,dwell, and retract
        
        Args:
            separation,zsnsr,force: offsets in their respective categories
        """
        self.approach.offset(separation,zsnsr,force)
        self.dwell.offset(separation,zsnsr,force)
        self.retract.offset(separation,zsnsr,force)
    def flip_forces(self):
        """
        multiplies all the forces by -1; useful after offsetting
        """
        self.approach.LowResData.force *= -1
        self.dwell.LowResData.force *= -1
        self.retract.LowResData.force *= -1
    def n_points_approach_dwell(self):
        """
        Returns:
            the number of points in the approach and dwell curves
        """
        return self.approach.Force.size + self.dwell.Force.size
    def get_retract_event_idx(self):
        """
        gets the slices of events *relative to the retract* (ie: idx 0 is
        the first point in the retract curve)
        
        Returns:
            list, each element is a slice like (start,stop,1) where start and   
            stop are the event indices
        """
        offset = self.n_points_approach_dwell() 
        # each event is a start/end tuple, so we just offset the min and max
        idx = [ slice(min(ev)-offset,max(ev)-offset,1) 
                for ev in self.retract.Events]
        return idx
    def get_retract_event_centers(self):
        """
        Returns:
            the mean of the event start and stop (its 'center')
        """
        get_mean = lambda ev: int(np.round(np.mean([ev.start,ev.stop]) ))
        return [ get_mean(ev) for ev in  self.get_retract_event_idx()]
    def get_predicted_retract_surface_index(self):
        """
        Assuming this have been zeroed, get the predicted retract surface index
        """
        filtered_obj = FEC_Util.GetFilteredForce(self.retract,
                                                 self.tau_num_points)
        return np.where(filtered_obj.Force >= 0)[0][0]

def get_surface_index(obj,n_smooth,last_less_than=True):
    """
    Get the surface index
    
    Args:
        obj: the timesepforce object to use
        n_smooth: number to smoothing   
        last_less_than: if true (default, 'raw' data), then we find the last
        time we are less than the baseline in obj.Force. Otherwise, the first
        time we are *greater* than...
    Returns 
        the surface index and baseline in force
    """
    force_baseline = np.median(obj.Force)
    filtered_obj = FEC_Util.GetFilteredForce(obj,n_smooth)
    if (last_less_than):
        # find the last time we are below the threshold ('raw' approach)
        search_func = lambda thresh: \
                np.where(filtered_obj.Force <= thresh)[0][-1]
    else:
        # find the first time we are above the threshhold ('processed' retract)
        search_func = lambda thresh: \
            np.where(filtered_obj.Force >= thresh)[0][0]
    idx_surface = search_func(force_baseline)
    # iterate once in order to get a better estimate of the baseline; we can
    # remove the effect of the invols entirely 
    if (last_less_than):
        force_baseline = np.median(obj.Force[:idx_surface])
    else:
        force_baseline = np.median(obj.Force[idx_surface:])
    idx_surface =  search_func(force_baseline)  
    return force_baseline,idx_surface,filtered_obj

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
    force_baseline,idx_surface,filtered_obj = \
        get_surface_index(approach,n_smooth,last_less_than=True)
    # get the separation at the baseline
    separation_baseline = filtered_obj.Separation[idx_surface]
    zsnsr_baseline = filtered_obj.Zsnsr[idx_surface]
    # zero everything 
    split_fec.zero_all(separation_baseline,zsnsr_baseline,force_baseline)
    if (flip_force):
        split_fec.flip_forces()
    split_fec.set_tau_num_points(n_smooth)
   
        
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

def spline_residual_mean_and_stdev(f,f_interp,start_q=1):
    """
    returns the mean and standard deviation associated with f-f_interp,
    from start_q% to 100-startq%
    
    Args:
        f: the 'noisy' function
        f_interp: the interpolated f (splined)
        start_q: the start perctile; we have to ignore huge outliers
    Returns:
        tuple of mean,standard deviation
    """
    # symetrically choose percentiles for the fit
    f_minus_mu = f-f_interp
    qr_1,qr_2 = np.percentile(a=f_minus_mu,q=[start_q,100-start_q])
    idx_fit = np.where( (f_minus_mu >= qr_1) &
                        (f_minus_mu <= qr_2))
    # fit a normal distribution to it, to get the standard deviation (globally)
    mu,std = norm.fit(f_minus_mu[idx_fit])
    return mu,std
    
def spline_gaussian_cdf(f,f_interp,std):
    """
    returns the CDF associated with the random variable with mean given by  
    f_interp and standard deviation associated with std, assuming gaussian
    about f-finterp
    
    Args:
        f: see spline_residual_mean_and_stdev
        f_interp: see spline_residual_mean_and_stdev
        std: standard deviation
    Returns:
        cummulative distribution
    """
    # get the distribution of the actual data
    distribution_force = norm(loc=f_interp, scale=std)
    # get the cdf of the data
    return distribution_force.cdf(f)
    
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

def auto_correlation_helper(auto):
    # normalize the auto correlation, add in a small bias to avoid 
    # taking the log of 0. data is normalized to 0->1, so it should be OK
    tol = 1e-9
    # auto norm goes from 0 to 1
    auto_norm = (auto - np.min(auto))/(np.max(auto)-np.min(auto)) 
    auto_median_normed = auto_norm - np.median(auto_norm)
    # statistical norm goes from -1 to 1
    statistical_norm = (auto_norm - 0.5) * 2
    log_norm = np.log(auto_norm + tol)
    fit_idx_max = np.where(auto_median_normed < 0.25)[0]
    assert fit_idx_max.size > 0 , "autocorrelation doesnt dip under threshold"
    # get the first time we cross under the threshold
    fit_idx_max =  fit_idx_max[0]
    return auto_norm,statistical_norm,log_norm,fit_idx_max
    
def auto_correlation_tau(x,f_user,deg_autocorrelation=1,fit_idx_max=None):
    """
    get the atucorrelation time of f (ie: fit polynomial to log(autocorrelation)
    vs x, so the tau is more or less the exponential decay constant
    
    Args:
        x: the unit of 'time'
        f_user: the function we want the autocorrelation of 
        deg_autocorrelation: the degree of autocorrelation to use. defaults to  
        linear, to get the 1/e time of autocorrelation
        
        fit_idx_max: maximum index to fit. defaults to until we hit 0 in
        the statistical autocorrelation 
    Returns:
        tuple of <autocorrelation tau, coefficients of log(auto) vs x fit,
                  auto correlation ('raw')>
    """
    f = f_user - np.mean(f_user)
    auto = np.correlate(f,f,mode='full')
    # only want the last half (should be identical?) 
    size = int(auto.size/2)
    auto = auto[size:]
    auto_norm,statistical_norm,log_norm,fit_idx_max_tmp = \
        auto_correlation_helper(auto)
    if fit_idx_max is None:
        fit_idx_max = fit_idx_max_tmp
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
    
