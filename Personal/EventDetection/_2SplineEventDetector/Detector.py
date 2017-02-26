# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import signal,stats

from Research.Personal.EventDetection.Util import Plotting
from GeneralUtil.python import PlotUtilities

def local_stdev(f,n):
    """
    Gets the local standard deviaiton (+/- n), except at boundaries 
    where it is just in the direction with data

    Args:
        f: what we want the stdev of
        n: window size (in either direction)
    Returns:
        array, same size as f, with the dat we want
    """
    max_n = f.size
    # go from (i-n to i+n)
    return np.array([np.std(f[max(0,i-n):min(max_n,i+n)]) 
                     for i in range(max_n)])

def adhesion_function_for_split_fec(split_fec):
    """
    returns the funciton adhesion_mask, with surface_index set to whatever
    the surface index of split_fec is predicted to be by the approach
    
    Args:
        split_fec: the split_force_extension object we want to mask the 
        adhesions of 
    Returns:
        lambda call to adhesion_mask, with arguments except for 
        probability_distribution,threshold supplied
    """
    surface_index = split_fec.get_predicted_retract_surface_index()
    n_points = split_fec.tau_num_points
    return (lambda *args,**kwargs: adhesion_mask(surface_index,n_points,
                                                 split_fec,
                                                 *args,**kwargs))
    
def _min_points_between(autocorrelation_tau_num_points):
    return int(np.ceil(autocorrelation_tau_num_points/2))
    
def adhesion_mask(surface_index,n_points,split_fec,
                  probability_distribution,threshold):
    """
    returns a boolean mask which is 0 where we can predict and adhesion 
    and zero elsewhere
    
    Args:
        surface_index: our best guess for where the surface is. 
        probability_distribution: see Detector._event_mask
        threshold: see Detector._event_mask
    Returns:
        list of event slices
    """
    to_ret = np.ones(probability_distribution.size,dtype=np.bool_)
    non_events = probability_distribution > threshold
    # determine the boundaries of the 'no events'
    min_points_between = _min_points_between(n_points)
    min_idx = surface_index + min_points_between    
    # remove all things before the predicted surface
    to_ret[:min_idx] = 0    
    no_event_mask = np.where(non_events)[0]
    # XXX finish current event, keep consuming events until startd/end
    # are beyond threshold
    event_mask = np.where(~non_events)[0]
    if (event_mask.size ==0 or no_event_mask.size == 0):
        return to_ret
    # POST: we have at least one event and one non-event 
    # (could be some adhesion!)
    # first, walk to where the smoothed y is at or abovethe median force
    # finally, make sure the smoothed force is back to zero
    # XXX switch to faster
    retract = split_fec.retract
    time = retract.Time
    smoothed_force = split_fec.retract_spline_interpolator()(time)
    # XXX fit a line to it instead?... by nice if this were approach based only
    force_threshold = np.median(smoothed_force[min_idx:])
    where_smoothed =  np.where(smoothed_force < force_threshold)[0]
    where_smoothed_and_greater = [e for e in where_smoothed if e >= min_idx] 
    if (len(where_smoothed_and_greater) == 0):
        return to_ret
    # POST: have some point where we are at or above the threshold
    min_idx = where_smoothed_and_greater[0]
    to_ret[:min_idx] = 0
    # determine events that contain the surface index
    event_boundaries = _event_slices_from_mask(event_mask,min_points_between)
    # get a list of the events with a starting point below the surface
    events_containing_surface = [e for e in event_boundaries  
                                 if (e.start <= min_idx)]
    if (len(events_containing_surface) == 0):
        return to_ret 
    # POST: at least one event contains the surface. Update the minimum index
    # to go to the end of the (last) event below or at the surface, unless
    # the end's end is below the surface, then just stick to our guns
    last_event_containing_surface_end = \
        events_containing_surface[-1].stop + min_points_between
    min_idx = max(min_idx,last_event_containing_surface_end)
    to_ret[:min_idx] = 0
    return to_ret                     
                     
class prediction_info:
    def __init__(self,event_idx,event_slices,local_stdev,interp,mask,
                 cdf,slice_fit,threshold,condition_result=None):
        """
        record the data from _predict_helper

        Args:
            event_idx : the event centers we found
            start_idx : the boundaries found for events
            local_stdev : the local standad deviation at each point in slice_fit
            mask : indices where an event is possible
            cdf : cummulative density function for the probability of a point
            in slice_fit, given the model of the data we have
        
            slice_fit : slice within the original retract where we tried to 
            find events. We have to remove the first few and last few points
            
            threshold: the threshhold for determining an event
            condition_result: mask used for adhesion, boolean array like cdf
        Returns:
            prediction_info object
        """
        self.event_idx = event_idx
        self.event_slices = event_slices
        self.local_stdev = local_stdev
        self.interp = interp
        self.mask = mask
        self.cdf = cdf
        self.slice_fit = slice_fit
        self.threshold = threshold
        self.condition_result = condition_result

def _event_mask(probability,threshold,condition_function=None):
    """
    Given a probability distribution and a threshold, returns the indices
    where the probability is less  than the threshhold
    
    Args;
        probability: array of numbers between 0 and 1
        threshold: maximum value in mask
        condition_function: a function taking in the probability and 
        threshold, and returning a boolean 1/0 array; a 1 is required for an 
        event
    Returns:
        tuple of (indices where the probability is less than the threshold
        and condition function is met, condition function boolean mask)
    """
    boolean_thresh = (probability <= threshold)
    if (condition_function is not None):  
        condition_result = condition_function(probability,threshold)
        conditions =(boolean_thresh & condition_result)
    else:
        condition_result = None
        conditions = boolean_thresh  
    return np.where(conditions)[0],condition_result
        
def _event_probabilities(x,y,interp,n_points,threshold):
    """
    determines the mask (and associated event detection information)
    
    Args:
        x,y: independent and dependent variable
        interp: the approximation to y vs x (ie: g*)
        n_points: number of points from autocorrelation function (ie: tau)

        threshold: maximum probability that a given datapoint fits the 
        model
    Returns:
        tuple of :
            probability_distribution : no-event probability for each point in y
            slice_fit : the part of x and y that mask is valid for
            stdevs: the local, windowed standard deviation, s(q)
    """
    min_points_between = _min_points_between(n_points)
    interp_first_deriv = interp.derivative(1)(x)
    # get the interpolated derivative
    interpolated_y = interp(x)
    # get a model for the local standard deviaiton using the autocorrelation
    # time from the event
    diff = y-interpolated_y
    stdevs = local_stdev(diff,n_points)
    # get the cwt of the wavelet; see pp219 of Mallat, Wavelet Tour (XXX TODO)
    median_local_stdev = np.median(stdevs)
    slice_fit = slice(min_points_between,-min_points_between,1)
    stdev_masked = stdevs[slice_fit]
    q25,q75 = np.percentile(stdev_masked,[25,75])
    iqr = q75-q25
    scale_idx = np.where( (stdev_masked <= q75) & (stdev_masked >= q25))
    scale = np.std(stdev_masked[scale_idx])
    # note: chebyshev is like
    # P(|X - mu| >=  k * sigma) <= 1/k^2
    # we write k = (s(q) - epsilon)/scale
    denom = (stdev_masked-median_local_stdev)
    k_chebyshev = denom/scale
    # note: chebyshev cant be more than 1 (could happen if the stdev is really 
    # close to the mean)
    chebyshev = np.minimum((1/k_chebyshev)**2,1)
    norm_dist = 1-stats.norm.cdf(stdev_masked,loc=median_local_stdev,
                                 scale=scale)
    # for the edge cases, assume the probability is one                         
    probability_distribution = np.ones(y.size)          
    # get the probability for all the non edge cases
    probability_distribution[slice_fit] = chebyshev
    return probability_distribution,slice_fit,stdevs

def _event_slices_from_mask(mask,min_points_between):
    """
    returns individual event slices for the mask, given that points shouldnt]
    be closer than min_points_between
    
    Args:
        mask: see _event_mask
        min_points_between: minimum number of points between events
        
    Returns:
        list of event slices
    """
    idx_step_changes = np.where(np.diff(mask) >= min_points_between)[0]
    # mask[idx_step_changes] gives the end of event i, for some i 
    # mask[idx_step_changes+1] gives the start of event (i+1), for the same i
    step_end_idx = mask[idx_step_changes]
    step_start_idx = mask[(idx_step_changes + 1)]
    # need to incllude the endpoints to get the proper events
    event_idx_end = list(step_end_idx) + [mask[-1]] 
    event_idx_start = [mask[0]] + list(step_start_idx)
    event_slices = [slice(start,end,1) 
                    for start,end in zip(event_idx_start,event_idx_end)]    
    return event_slices

def _loading_rate_helper(x,y,slice_event):
    """
    Determine where a (single, local) event is occuring in the slice_event
    (of length N) part of x,y by:
    (1) Finding the maximum of y in the slice
    (2) Fitting a line to the N points up to the maximum
    (3) Determining the last point at which y[slice_event] is above the 
    predicted line from (2). If this doesnt exist, just uses the maximum

    Args:
        x, y: x and y values. we assume an event is from high to low in y
        slice_event: where to fit
    Returns:
        tuple of <fit_x,fit_y,predicted y based on fit, idx_above_predicted>
    """
    # determine the local maximum
    offset = slice_event.start
    n_points = int(np.ceil((slice_event.stop-offset+1)/2))
    y_event = y[slice_event]
    x_event = x[slice_event]
    local_max_idx = offset + np.argmax(y_event)
    n = x.size
    # end our fit at the midpoint of the event; start offset from 
    end_fit_idx = local_max_idx
    start_fit_idx  = max(0,end_fit_idx-n_points)
    fit_slice = slice(start_fit_idx,end_fit_idx,1)
    # fit 1-D until the local max
    fit_x = x[fit_slice]
    fit_y = y[fit_slice]
    coeffs = np.polyfit(x=fit_x,y=fit_y,deg=1)
    pred = np.polyval(coeffs,x=x_event)
    # determine where the data *in the __original__ slice* is __last__
    # above the fit (after that, it is consistently below it)
    idx_above_predicted_rel = np.where(y_event > pred)[0]
    # dont look at things past where we fit...
    idx_above_predicted = [offset + i for i in idx_above_predicted_rel]
    return fit_x,fit_y,pred,idx_above_predicted
    
def event_by_loading_rate(*args,**kwargs):
    """
    see _loading_rate_helper 

    Args:
        see _loading_rate_helper
    Returns:
        predicted index (absolute) in x,y where we think the event is happening
    """
    fit_x,fit_y,pred,idx_above_predicted = _loading_rate_helper(*args,**kwargs)
    # POST: have a proper max, return the last time we are above
    # the linear prediction
    if (len(idx_above_predicted) == 0):
        return fit_max_idx
    return idx_above_predicted[-1]

def _predict(x,y,n_points,interp,threshold,local_event_idx_function,
             condition_function=None):
    """
    general method to predict the event boundaries and centers
    
    Args:
        x: see _event_probabilities
        y: see _event_probabilities
        n_points: see _event_probabilities
        interp: see _event_probabilities
        threshold: see _event_probabilities
        local_event_idx_function: a function which takes a slice of x,y,slice
        as its  only argument and returns the most likely index of an event. 
        the slice passsed should have only one event
        
        condition_function: see _event_mask
    Returns:
        list of event slices
    """
    min_points_between = _min_points_between(n_points)
    probability_distribution,slice_fit,stdevs = \
        _event_probabilities(x,y,interp,n_points,threshold)
    mask,condition_result = _event_mask(probability_distribution,
                                        threshold,condition_function)
    # only keep points where we are farther than min_points between from the 
    # edges (ie: from index 0 and N-1)
    n = mask.size
    mask = mask[np.where( (mask >= min_points_between ) | 
                          (mask <= n-min_points_between ) )]
    if (mask.size > 0):
        event_slices = _event_slices_from_mask(mask,min_points_between)
    else:
        event_slices = []
    # determine where the events are happening locally (guarentee at least
    # a search window of min_points)
    remainder_split = [ int(np.ceil((min_points_between-(e.stop-e.start)/2)))
                        for e in event_slices]
    event_slices = [slice(event.start-remainder,event.stop+remainder,1) 
                    for event,remainder in zip(event_slices,remainder_split)]
    event_idx = [local_event_idx_function(x,y,e) for e in event_slices]
    to_ret = prediction_info(event_idx = event_idx,
                             event_slices = event_slices,
                             local_stdev = stdevs,
                             interp = interp,
                             mask = mask,
                             cdf=probability_distribution,
                             slice_fit=slice_fit,
                             threshold=threshold,
                             condition_result=condition_result)
    return to_ret                                
                             
def _predict_helper(split_fec,threshold,**kwargs):
    """
    uses spline interpolation and local stadard deviations to predict
    events.

    Args:
        split_fec: split_force_extension object, already initialized, and 
        zerod, with the autocorraltion time set. 

        threshhold: maximum probability that a given datapoint fits the 
        model
        
        kwargs: passed to _predict
    Returns:
        prediction_info object
    """
    retract = split_fec.retract
    time,separation,force = retract.Time,retract.Separation,retract.Force
    n_points = split_fec.tau_num_points
    local_event_idx_function = \
        lambda x,y,event : event.start + np.argmax(y[event])
    # N degree b-spline has continuous (N-1) derivative
    interp = split_fec.retract_spline_interpolator(deg=2)
    to_ret = _predict(x=time,
                      y=force,
                      n_points=n_points,
                      interp=interp,
                      threshold=threshold,
                      local_event_idx_function=event_by_loading_rate,
                      **kwargs)
    # XXX modify mask; find first time under threshhold after where we predict
    # the surface
    return to_ret
