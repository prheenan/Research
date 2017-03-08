# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import signal,stats

from Research.Personal.EventDetection.Util import Analysis,Plotting
from GeneralUtil.python import PlotUtilities,GenUtilities

from scipy.ndimage.filters import uniform_filter1d



from itertools import chain

def join_contiguous_slices(slices, offset=0):
    """
    # XXX see:
    stackoverflow.com/questions/24317211/
    merge-overlapping-numeric-ranges-into-continuous-ranges
    """
    flatten = chain.from_iterable
    LEFT, RIGHT = 1, -1
    data = [[s.start,s.stop] for s in slices]
    data = sorted(flatten(((start, LEFT), (stop + offset, RIGHT))
            for start, stop in data))
    c = 0
    for value, label in data:
        if c == 0:
            x = value
        c += label
        if c == 0:
            yield slice(x, value - offset,1)

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
    """
    for linear stdev, see: 
    stackoverflow.com/questions/18419871/
    improving-code-efficiency-standard-deviation-on-sliding-windows
    """
    mode = 'reflect'
    c1 = uniform_filter1d(f, size=n*2, mode=mode, origin=-n)
    c2 = uniform_filter1d(f*f, size=n*2, mode=mode, origin=-n)
    # sigma^2 = ( <x^2> - <x>^2 )^(1/2), shouldnt dip below 0
    safe_variance = np.maximum(0,c2 - c1*c1)
    stdev = (safe_variance**.5)
    return stdev

def spline_derivative_probability(split_fec):
    """
    see mask_spline_derivative, except this returns the probability 
    distribution of each point p, P(deriv to be less than p) by chebyshev
    
    Args:
        split_fec: the split_force_extension object we want to mask the 
        adhesions of 
    Returns:
        array, 1 where we are within a stdev of the median, otherwise
        (stdev/(p-std))**2
    """
    retract = split_fec.retract
    time = retract.Time 
    interpolator = split_fec.retract_spline_interpolator()
    return _spline_derivative_probability_generic(time,interpolator)

def _spline_derivative_probability_generic(x,interpolator):
    """
    see  spline_derivative_probability, except a genertic method
    
    Args:
        x: x values
        interpolator: to interpolate along
    Returns:
        see spline_derivative_probability
    """
    derivative_force = interpolator.derivative()(x)
    # get the median and std of deriv
    med_deriv = np.median(derivative_force)
    q_loq_percentile = 0
    q_high_percentile = 100
    q_low,q_high = np.percentile(derivative_force,
                                 [q_loq_percentile,q_high_percentile])
    iqr_region_idx = np.where( (derivative_force <= q_high) & 
                               (derivative_force >= q_low))[0]
    std_iqr = np.std(derivative_force[iqr_region_idx])
    probability = np.zeros(derivative_force.size)
    # anything at or above the median isnt interesting
    probability[np.where(derivative_force >= med_deriv - std_iqr)]  = 1
    # other things might be
    possible_idx = np.where(derivative_force < med_deriv)
    possible_deriv = derivative_force[possible_idx]
    k = (possible_deriv-med_deriv)/std_iqr
    probability[possible_idx]  = 1/k**2
    probability = np.minimum(probability,1)
    return probability 

def mask_spline_derivative(split_fec):
    """
    returns a mask on the spline derivative: we must be at least one 
    standard deviation away for the mask to be one
    
    Args:
        split_fec: the split_force_extension object we want to mask the 
        derivative of 
    Returns:
        0/1 array of the same size as split_fec.retract.Force
    """
    return (spline_derivative_probability(split_fec) < 1)

def derivative_mask_function(split_fec,slice_to_use,
                             boolean_array,probability,threshold,
                             *args,**kwargs):
    """
    returns mask_spline_derivative
    
    Args:
        split_fec: the split_force_extension object we want to mask the 
        derivative of 

        *args,**kwargs: ignored
    Returns:
        see adhesion_mask_function_for_split_fec, except derivative mask
    """
    slice_v = slice_to_use
    offset = slice_v.start    
    n_points = split_fec.tau_num_points
    min_points_between = _min_points_between(n_points)    
    # if we dont have any points, return
    if ((slice_v.stop - offset) < min_points_between):
        return slice_to_use,boolean_array,probability
    # POST: something to look at. find the spline-interpolated derivative
    # probability
    retract = split_fec.retract
    x = retract.Time[slice_v]
    interp = split_fec.retract_spline_interpolator(slice_to_fit=slice_v)
    spline_probability_in_slice=\
        _spline_derivative_probability_generic(x,interp)
    # determine where the derivative is possibly outlying; that is a necessary
    # but not sufficient condition for an event
    tol = 1e-9
    possible_event =  (1 - spline_probability_in_slice) > tol
    event_slice_mask = np.where(possible_event)[0]
    if (event_slice_mask.size == 0):
        return slice_to_use,boolean_array,probability
    # POST: have at least one possible spline event
    # only consider events which start after the offset, to avoid edge effects
    absolute_min_idx = offset + min_points_between          
    event_boundaries = _event_slices_from_mask(event_slice_mask,
                                               min_points_between)    
    events_starting_before_min_idx = [e for e in event_boundaries 
                                      if e.start < min_points_between]
    if (len(events_starting_before_min_idx) > 0):
        event_end_in_slice = events_starting_before_min_idx[0].stop
        # XXX Wtf?
        absolute_min_idx = max(absolute_min_idx,event_end_in_slice+offset)
    # POST: absolute_min_idx is the index *in the original array* where we 
    # should start looking for events. 
    probability_updated = probability.copy()
    # determine where the probability is
    where_no_event = np.where(~possible_event)
    spline_boolean = np.zeros_like(boolean_array)
    spline_boolean[slice_v] = possible_event
    probability_updated[slice_v] *= spline_probability_in_slice
    probability_updated[slice_v][where_no_event] = 1
    # mask everything until the minimum index; probability is set to one,
    # boolean mask is set to zero
    spline_boolean[:absolute_min_idx] = 0
    probability_updated[:absolute_min_idx] = 1
    slice_updated = slice(absolute_min_idx,slice_to_use.stop,1)
    return slice_updated,(spline_boolean & boolean_array),probability_updated

def adhesion_mask_function_for_split_fec(split_fec,slice_to_use,boolean_array,
                                         probability,threshold,
                                         *args,**kwargs):
    """
    returns the funciton adhesion_mask, with surface_index set to whatever
    the surface index of split_fec is predicted to be by the approach
    
    Args:
        split_fec: the split_force_extension object we want to mask the 
        adhesions of 

        *args,**kwargs: see adhesion
    Returns:
        new mask and probability distribution
    """
    surface_index = split_fec.get_predicted_retract_surface_index()
    n_points = split_fec.tau_num_points
    bool_adhesion = adhesion_mask(surface_index,n_points,split_fec,
                                  probability_distribution=probability,
                                  threshold=threshold,*args,**kwargs)
    where_not_adhesion = np.where(bool_adhesion)[0]
    if (where_not_adhesion.size > 1):
        adhesions_removed = boolean_array.copy()
        start = where_not_adhesion[0]
        end = where_not_adhesion[-1]
        # remove the starting and ending adhesion
        adhesions_removed[:start] = 0
        adhesions_removed[end:] = 0
        # update the no-event probabilities
        slice_update = slice(start,end,1)
        retract = split_fec.retract
        time = retract.Time
        force = retract.Force
        interp = split_fec.\
            retract_spline_interpolator(slice_to_fit=slice_update)
        probability_updated = np.ones( boolean_array.size)
        prob_tmp, _ = _no_event_probability(time,interp,force,n_points,
                                            slice_fit=slice_update)
        probability_updated[slice_update] = np.minimum(1,prob_tmp)
    else:
        slice_update = slice(0,None,1)
    return slice_update,boolean_array,probability

def _min_points_between(autocorrelation_tau_num_points):
    """
    returns the minimum recquired points between two discrete events,
    given a number of filtering points
    
    Args:
        autocorrelation_tau_num_points: number of filtering points
    """
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
    # remove all things before the predicted surface, and at the boundary
    to_ret[:min_idx] = 0    
    to_ret[-min_idx:] = 0    
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
    # the end is below the surface, then just keep the minimum at the surface
    last_event_containing_surface_end = \
        events_containing_surface[-1].stop + min_points_between
    min_idx = max(min_idx,last_event_containing_surface_end)
    to_ret[:min_idx] = 0
    return to_ret                     
                     

class prediction_info:
    def __init__(self,event_idx,event_slices,local_stdev,interp,mask,
                 cdf,slice_fit,threshold,probabilities=None,
                 condition_results=None):
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
            condition_results: list of masks used for adhesion, 
            list of boolean arrays, one per mask 
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
        self.condition_results = condition_results
        self.probabilities= probabilities

def _event_mask(probability,threshold):
    """
    Given a probability distribution and a threshold, returns the indices
    where the probability is less  than the threshhold
    
    Args:
        example_split: the split force extension curve to use
        probability: array of numbers between 0 and 1
        threshold: maximum value in mask
        condition_functions: list of functions taking in the probability and 
        threshold, and returning a boolean 1/0 array; a 1 is required for an 
        event
    Returns:
        tuple of (indices where the probability is less than the threshold)
    """
    boolean_thresh = (probability <= threshold)
    return np.where(boolean_thresh)[0]

def _no_event_probability(x,interp,y,n_points,slice_fit=None):
    """
    returns the no-event probability at each point in y

    Args:
        x: the x values that interp takes, see _event_probabilities
        y: the y values we are searching for an event, see _event_probabilities
        interp: see _event_probabilities
        n_points: number of points to use in estimating r(q)=g-g* by the 
        local standard deviaiton of y-interp(x)
    
        slice_fit: an optional slice to use to compute the probabilities
    Returns:
        tuple of <probability, local stdevs>
    """
    if (slice_fit is None):
        slice_fit = slice(0,None,1)
    x = x[slice_fit]
    y = y[slice_fit]
    # get the interpolated function
    interpolated_y = interp(x)
    # get a model for the local standard deviaiton
    diff = y-interpolated_y
    stdevs = local_stdev(diff,n_points)
    # get the cwt of the wavelet; see pp219 of Mallat, Wavelet Tour (XXX TODO)
    median_local_stdev = np.median(stdevs)
    stdev_masked = stdevs[slice_fit]
    qlow,qhigh = np.percentile(stdev_masked,[0,75])
    iqr = qlow-qhigh
    scale_idx = np.where( (stdev_masked >= qlow) & (stdev_masked <= qhigh))
    scale = np.std(stdev_masked[scale_idx])
    # note: chebyshev is like
    # P(|X - mu| >=  k * sigma) <= 1/k^2
    # we write k = (s(q) - epsilon)/scale
    denom = (stdev_masked-median_local_stdev)
    k_chebyshev = denom/scale
    # determine where the chebyshev is 'safe', otherwise we are at or above
    # the mean estimate and hence not a useful metric
    cheby_idx = np.where(k_chebyshev >= 1)
    chebyshev = np.ones(k_chebyshev.size)
    # actually calculate the upper bound for the probability
    chebyshev[cheby_idx] = (1/k_chebyshev[cheby_idx])**2
    # for the edge cases, assume the probability is one                         
    probability_distribution = np.ones(y.size)      
    # get the probability for all the non edge cases
    probability_distribution[slice_fit] = chebyshev
    return probability_distribution,stdevs
        
def _event_probabilities(x,y,interp,n_points,threshold):
    """
    determines the mask (and associated event detection information)
    
    Args:
        x,y: independent and dependent variable (ie: 'q' and 'g'
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
    slice_fit = slice(min_points_between,-min_points_between,1)
    probability_distribution_in_slice,stdevs = \
        _no_event_probability(x,interp,y,n_points=n_points,slice_fit=slice_fit)
    probability_distribution = np.ones_like(x)
    probability_distribution[slice_fit] =  probability_distribution_in_slice
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
    return fit_x,fit_y,pred,idx_above_predicted,local_max_idx
    
def event_by_loading_rate(*args,**kwargs):
    """
    see _loading_rate_helper 

    Args:
        see _loading_rate_helper
    Returns:
        predicted index (absolute) in x,y where we think the event is happening
    """
    fit_x,fit_y,pred,idx_above_predicted,local_max_idx = \
            _loading_rate_helper(*args,**kwargs)
    # POST: have a proper max, return the last time we are above
    # the linear prediction
    if (len(idx_above_predicted) == 0):
        return local_max_idx
    return idx_above_predicted[-1]


def _predict(x,y,n_points,interp,threshold,local_event_idx_function,
             remasking_functions=None):
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
        
        remasking_functions: XXX
    Returns:
        list of event slices
    """
    min_points_between = _min_points_between(n_points)
    probability_distribution,slice_fit,stdevs = \
        _event_probabilities(x,y,interp,n_points,threshold)
    bool_array = probability_distribution < threshold
    masks = [np.where(bool_array)[0]]
    probabilities = [probability_distribution.copy()]
    slice_to_use = slice(0,None,1)
    if (remasking_functions is not None):
        for f in remasking_functions:
            res = f(slice_to_use=slice_to_use,
                    boolean_array=bool_array,
                    probability=probability_distribution,
                    threshold=threshold)
            slice_to_use,bool_array, probability_distribution = res
            probabilities.append(probability_distribution)
            masks.append(np.where(bool_array)[0])
    # only keep points where we are farther than min_points between from the 
    # edges (ie: from index 0 and N-1)
    mask = np.where(bool_array)[0]
    n = mask.size
    mask = mask[np.where( (mask >= min_points_between ) | 
                          (mask <= n-min_points_between ) )]
    if (mask.size > 0):
        event_slices = _event_slices_from_mask(mask,min_points_between)
    else:
        event_slices = []
    # XXX reject events with a very small time?
    event_duration = [ (e.stop-e.start) for e in event_slices]
    # determine where the events are happening locally (guarentee at least
    # a search window of min_points)
    remainder_split = [ int(np.ceil((min_points_between-(e.stop-e.start)/2)))
                        for e in event_slices]
    event_slices = [slice(event.start-remainder,event.stop+remainder,1) 
                    for event,remainder in zip(event_slices,remainder_split)]
    # POST: slices are of length min points, determine which events overlap, 
    # combine them if they do.
    event_slices = list(join_contiguous_slices(event_slices))
    # POST: event slices aren't contiguous
    event_idx = [local_event_idx_function(x,y,e) for e in event_slices]
    to_ret = prediction_info(event_idx = event_idx,
                             event_slices = event_slices,
                             local_stdev = stdevs,
                             interp = interp,
                             mask = mask,
                             cdf=probability_distribution,
                             slice_fit=slice_fit,
                             threshold=threshold,
                             condition_results=masks,
                             probabilities=probabilities)
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

def _predict_functor(example,f):
    """
    python doesn't like creating lambda functions using 
    function references in a list (2017-3-1), so I use a functor instead

    Args:
        example: first argument of type like f. split_fec is used in predict
        f: function we are calling. should take split_fec, then *args,**kwargs
        (see _predict
    returns:
        a lambda function passing arguments and keyword argument to f 
    """
    return lambda *args,**kwargs : f(example,*args,**kwargs)


def _predict_full(example,threshold=1e-2):
    """
    see predict, example returns tuple of <split FEC,prediction_info>
    """
    example_split = Analysis.zero_and_split_force_extension_curve(example)
    f_refs = [adhesion_mask_function_for_split_fec,derivative_mask_function]
    funcs = [ _predict_functor(example_split,f) for f in f_refs]
    final_dict = dict(remasking_functions=funcs,
                      threshold=threshold)
    pred_info = _predict_helper(example_split,**final_dict)
    return example_split,pred_info

def predict(example,threshold=1e-2):
    """
    predict a single event from a force extension curve

    Args:
        example: TimeSepForce
        threshold: maximum probability under the no-event hypothesis
    Returns:
        list of event starts
    """
    example_split,pred_info = _predict_full(example,threshold=threshold)
    return pred_info.event_idx

