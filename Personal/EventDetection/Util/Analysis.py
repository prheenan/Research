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
    def retract_spline_interpolator(self,slice_to_fit=None,**kwargs):
        """
        returns an interpolator for force based on the stored time constant tau
        for the retract force versus time curbe

        Args:
            slice_to_fit: which part of the retract to fit
            kwargs: passed to spline_interpolator
        """
        return spline_fit_fec(self.tau,
                              self.retract,slice_to_fit=slice_to_fit,**kwargs)
    def approach_spline_interpolator(self,slice_to_fit=None,**kwargs):
        """
        See retract_spline_interpolator, but for the approach
        """
        return spline_fit_fec(self.tau,self.approach,
                              slice_to_fit=slice_to_fit,**kwargs)        
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
    def surface_distance_from_trigger(self):
        """
        returns the distance in separtion units from the trigger point
        """
        return abs(min(self.approach.Separation))
    def get_predicted_approach_surface_index(self):
        """
        returns the predicted place the surface is on the approach
        """    
        return np.where(self.approach.Force >0)[0][-1]
    def get_predicted_retract_surface_index(self):
        """
        Assuming this have been zeroed, get the predicted retract surface index
        """
        n_points = _index_surface_relative(self.retract.Separation,
                                           self.surface_distance_from_trigger())
        return n_points

def _index_surface_relative(x,offset_needed):
    """
     returns a crude estimate for  the predicted index offset for the surface
        
    Args:
        x: the time series of separation
        offset_needed: the x offset 
    Returns: 
        number of points for x to displace by offset_needed
    """    
    sep_diff = np.median(np.abs(np.diff(x)))
    n_points = int(np.ceil(offset_needed/sep_diff))
    return n_points
        
def spline_fit_fec(tau,time_sep_force,slice_to_fit=None,**kwargs):
    """
    returns an interpolator object on the given TimeSepForce object
     
    Args:
        tau: see spline_interpolator
        time_sep_force: get t he time and force from this as x,y to 
        spline_interpolator
        
        slice_to_fit: part of x,f to fit
        **kwargs: passed to spline_interpolator
    returns:
        see spline_interpolator
    """    
    x,f = time_sep_force.Time,time_sep_force.Force
    if (slice_to_fit is None):
        slice_to_fit = slice(0,None,1)
    return spline_interpolator(tau,x[slice_to_fit],f[slice_to_fit],
                               **kwargs)        
        
def filter_fec(obj,n_points):
    return FEC_Util.GetFilteredForce(obj,n_points,spline_interpolated_by_index)

def bhattacharyya_probability_coefficient_1d(v1,v2,bins):
    """
    # return the bhattacharyya distance between two 1-d arras

    Args:
        v<1/2>: see  bhattacharyya_probability_coefficient_dd, except 1-D      
        bins: how to bin them, see 
    Returns:
        bhattacharyya distance, see bhattacharyya_probability_coefficient
    """
    return bhattacharyya_probability_coefficient_dd(v1,v2,[bins])

def bhattacharyya_probability_coefficient_dd(v1,v2,bins):
    """
    # return the bhattacharyya distance between arbitrary-dimensional
    #probabilities, see  bhattacharyya_probability_coefficient

    Args:
        v<1/2>: two arbitrary-dimensional lists to compare
        bins: how to bin them
    Returns:
        bhattacharyya distance, see bhattacharyya_probability_coefficient
    """
    histogram_kwargs = dict(bins=bins,weights=None,normed=False)
    v1_hist,v1_edges = np.histogramdd(sample=v1,**histogram_kwargs)
    v2_hist,v2_edges = np.histogramdd(sample=v2,**histogram_kwargs)
    return bhattacharyya_probability_coefficient(v1_hist,v2_hist)

def bhattacharyya_probability_coefficient(v1_hist,v2_hist):
    """
    # return the bhattacharyya distance between the probabilities, see:
    # https://en.wikipedia.org/wiki/Bhattacharyya_distance

    Args:
        v<1/2>_hist: values of two ditributions in each bins
    Returns:
        bhattacharyya distance
    """
    v1_hist = v1_hist.flatten()
    v2_hist = v2_hist.flatten()
    p1 = v1_hist/np.sum(v1_hist)
    p2 = v2_hist/np.sum(v2_hist)
    return sum(np.sqrt(p1*p2))
    

def _surface_index(filtered_y,y,last_less_than=True):
    """
    Get the surface index
    
    Args:
        y: the y we are searching for the surface of (raw)
        filtered_y: the filtered y value
        n_smooth: number to smoothing   
        last_less_than: if true (default, 'raw' data), then we find the last
        time we are less than the baseline in obj.Force. Otherwise, the first
        time we are *greater* than...
    Returns 
        the surface index and baseline in force
    """
    median = np.median(y)
    lt = np.where(y < median)[0]
    # determine the last time we were less than the median;
    # use this as a marker between the invols and the surface region
    last_lt = lt[-1]
    x = np.arange(start=0,stop=y.size,step=1)
    x_approach = x[:last_lt]
    x_invols = x[last_lt:]
    coeffs_approach = np.polyfit(x=x_approach,y=y[:last_lt],deg=1)
    coeffs_invols = np.polyfit(x=x_invols,y=y[last_lt:],deg=1)
    pred_approach = np.polyval(coeffs_approach,x=x)
    pred_invols = np.polyval(coeffs_invols,x=x)
    surface_idx = np.argmin(np.abs(pred_approach-pred_invols))
    return median,surface_idx

def get_surface_index(obj,n_smooth,last_less_than=True):
    """
    Get the surface index
    
    Args:
        see _surface_index
    Returns 
        see _surface_index, except last (extra) tuple element is filtered
        obj 
    """
    filtered_obj = filter_fec(obj,n_smooth)
    baseline,idx = _surface_index(filtered_obj.Force,obj.Force,
                                  last_less_than=last_less_than)
    return baseline,idx,filtered_obj

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
    
def spline_interpolated_by_index(f,nSmooth,**kwargs):
    """
    returnsa spline interpolator of f versus 0,1,2,...,(N-1)
    
    Args:
        f: function to interpolate
        nSmooth: distance between knots (smoothing number)
        **kwargs: passed to spline_interpolator
    Returns: 
        spline interpolated value of f on the indices (*not* an interpolator
        object, just an array) 
    """
    x = np.arange(start=0,stop=f.size,step=1)
    return spline_interpolator(nSmooth,x,f,**kwargs)(x)

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
    
def auto_correlation_tau(x,f_user,deg_autocorrelation=1,
                         autocorrelation_skip_points=None,fit_idx_max=None):
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
    if (autocorrelation_skip_points is not None):
        auto = auto[autocorrelation_skip_points:]
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

def zero_and_split_force_extension_curve(example):
    """
    zeros a force extension curve by its meta information and the touchoff
    on the approach

    Args:
        example: 'raw' force extension to use (negative force is away
        from surface on molecule)
    returns:
        example as an Analysis.split_force_extension object
    """
    example_split = split_FEC_by_meta(example)
    approach = example_split.approach
    retract = example_split.retract 
    f = approach.Force
    x = approach.Time
    num_points = int(np.ceil(f.size * 0.025))
    # zero out everything to the approach using the autocorrelation time 
    zero_by_approach(example_split,num_points)
    return example_split

def loading_rate_rupture_force_and_index(time,force,slice_to_fit):
    """
    given a portion of time and force to fit, the loading rate is determined 
    by the local slope. The rupture force is determined by finding the last
    time the (XXX should use two lines in case of flickering)
    
    Args:
        time/force: should be self-explanatory. Force should be zeroed.
        slice_to_fit: where we are fitting
    Returns:
        tuple of <loading rate,rupture force,index_of_rupture>
    """
    x = time[slice_to_fit]
    y = force[slice_to_fit]
    # XXX can fit a line, throw an error?
    if (x.size < 2):
        raise IndexError("Can't fit a line to something with <2 points")
    coeffs = np.polyfit(y=y,x=x,deg=1)
    predicted = np.polyval(coeffs,x=x)
    loading_rate, _ = coeffs
    # determine the last time the *data* is above the prediction
    where_above = np.where(y > predicted)[0]
    if (where_above.size == 0):
        # unlikely but worth checking
        last_idx_above = np.argmax(y)
    else:
        last_idx_above = where_above[-1]
    # determine what we *predict* to be the value at that point
    rupture_force = predicted[last_idx_above]
    return loading_rate,rupture_force,last_idx_above
    
    
