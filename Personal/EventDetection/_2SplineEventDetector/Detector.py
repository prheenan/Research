# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import signal,stats

from Research.Personal.EventDetection.Util import Analysis,Plotting
from GeneralUtil.python import PlotUtilities,GenUtilities
# XXX reduce import size below
from Research.Personal.EventDetection._2SplineEventDetector._no_event import \
    _min_points_between,_predict,\
    _probability_by_cheby_k,_no_event_chebyshev,_event_slices_from_mask
from Research.Personal.EventDetection._2SplineEventDetector import _no_event

def get_slice_by_max_value(interp_sliced,offset,slice_list):
    value_max = [max(interp_sliced[e.start-offset:e.stop-offset])
                 for e in slice_list]
    return np.argmax(value_max)

def force_value_mask_function(split_fec,slice_to_use,
                              boolean_array,probability,threshold,
                              *args,**kwargs):
    """
    masks the interpolated force to be at least one r(q) above the median

    Args:
         see adhesion_mask_function_for_split_fec 
    Returns:
         see adhesion_mask_function_for_split_fec, except tuples are dealing
         with the force value mask.
    """
    retract = split_fec.retract
    n_points = split_fec.tau_num_points
    min_points_between = _min_points_between(n_points)
    f = retract.Force[slice_to_use]
    x = retract.Time[slice_to_use]
    boolean_ret = boolean_array.copy()
    probability_updated = probability.copy()
    interpolator = split_fec.retract_spline_interpolator(slice_to_use)
    interp_f = interpolator(x)
    diff = f-interp_f
    stdev = Analysis.local_stdev(diff,n_points)
    med = np.median(interp_f)
    epsilon,sigma = split_fec.get_epsilon_and_sigma()
    # essentially: when is the interpolated value 
    # at least one (local) standard deviation above the median
    # we admit an event might be possible
    thresh = sigma
    local_integral = Analysis.local_integral(stdev-epsilon,min_points_between)
    # the threshold is the noise sigma times  the number of points 
    # (2*num_between)
    thresh_integral = 2 * sigma * min_points_between
    probability_updated[slice_to_use] *= \
            _no_event_chebyshev(local_integral,0,thresh_integral)
    boolean_ret[slice_to_use] *= (probability_updated[slice_to_use] < threshold)
    boolean_updated,probability_updated = ret
    return slice_to_use,boolean_updated,probability_updated

def safe_reslice(original_boolean,original_probability,condition,
                 min_points_between,get_best_slice_func):
    """
    applies the boolean <condition> array without creating new events; if 
    more than one event exists where a previous one existed, picks by 
    get_best_slice_func

    Args;
        original_boolean: boolean array, size N
        original_probability: probability array, size N
        condition: array to apply to above, size N. where condition is 1,
        the boolean array becomes *zero* 

        min_points_between: see _event_slices_from_mask
        get_best_slice_func: given a list of slices in the *original* data
        (ie: size N), this function should return the index of the single
        slice to keep
    Returns:
        tuple of <updated boolean, updated probability>
    """
    new_boolean = original_boolean.copy()
    new_probability = original_probability.copy()
    where_condition = np.where(condition)[0]
    if (where_condition.size > 0):
        new_boolean[where_condition] = 0 
        new_probability[where_condition] = 1
        # get the original and new slices
        mask_original = np.where(original_boolean)[0]
        mask_new = np.where(new_boolean)[0]
        if (mask_new.size == 0 or mask_original.size == 0):
            return new_boolean,new_probability
        # POST: have something to do
        original_events = _event_slices_from_mask(mask_original,
                                                  min_points_between)
        new_events = _event_slices_from_mask(mask_new,min_points_between)
        remove,keep = [],[]
        for e in original_events:
            start,stop = e.start,e.stop
            # determine which new events are within the old events
            candidates = [tmp for tmp in new_events 
                          if (tmp.start >= start) and (tmp.stop <= stop)
                          and (tmp.start < tmp.stop)]
            if (len(candidates) == 0):
                continue
            # determine the best new event within the old event, in the subslice
            # indices
            idx_best = get_best_slice_func(candidates)
            # get the best
            keep.append(candidates[idx_best])
            # remove the events
            remove.extend(candidates[i] for i in range(len(candidates))
                          if i != idx_best)
        # anything left over should also be deleted
        remove.extend(new_events)
        # POST: know what to remove and what to keep
        new_boolean = np.zeros(new_boolean.size)
        new_probability = np.ones(new_probability.size)
        for keep_idx in keep:
            new_boolean[keep_idx] = 1
            new_probability[keep_idx] = original_probability[keep_idx]
        # pick out the minimum derivative slice within each previous slice
    return new_boolean,new_probability

def derivative_mask_function(split_fec,slice_to_use,
                             boolean_array,probability,threshold,
                             no_event_parameters_object,
                             negative_only=True,
                             *args,**kwargs):
    """
    returns : see adhesion_mask_function_for_split_fec 
    
    Args:
        split_fec: the split_force_extension object we want to mask the 
        derivative of 
    
        other arguments: see adhesion_mask
        *args,**kwargs: ignored
    Returns:
        see adhesion_mask_function_for_split_fec, except derivative mask
    """
    offset = slice_to_use.start    
    n_points = split_fec.tau_num_points
    min_points_between = _min_points_between(n_points)    
    # if we dont have any points, return
    boolean_ret = boolean_array.copy()
    probability_updated = probability.copy()
    n_full_force = split_fec.retract.Force.size
    n_full_force = split_fec.retract.Force.size
    stop = n_full_force if (slice_to_use.stop) is None else slice_to_use.stop
    if ((stop - offset) < min_points_between):
        return slice_to_use,boolean_ret,probability_updated
    # POST: something to look at. find the spline-interpolated derivative
    # probability
    retract = split_fec.retract
    time = retract.Time
    force = retract.Force
    x_sliced =  time[slice_to_use]
    interp = split_fec.retract_spline_interpolator(slice_to_fit=slice_to_use)
    interp_deriv = interp.derivative()(x_sliced) 
    # POST: start looking at other points
    median = np.median(interp_deriv)
    # get rid of final outlying derivative points 
    where_above = np.where(interp_deriv < median)[0]
    where_below = np.where(interp_deriv > median)[0]
    last_index = offset + min(where_above[-1],where_below[-1])
    absolute_min_idx = offset +min_points_between
    absolute_max_index = min(slice_to_use.stop,last_index)
    # remove the bad points
    boolean_ret[:absolute_min_idx] = 0
    boolean_ret[absolute_max_index:] = 0
    slice_to_use = slice(absolute_min_idx,absolute_max_index,1)
    x_sliced =  time[slice_to_use]
    # get the derivative probability
    probability_deriv = \
        _no_event._derivative_probability(interp,x_sliced,
                                          no_event_parameters_object)
    # modulate the probabilities by the approach
    probability_updated[slice_to_use] *= probability_deriv
    boolean_ret[slice_to_use] = probability_updated[slice_to_use] < threshold
    return slice_to_use,boolean_ret,probability


def integral_mask_function(split_fec,slice_to_use,
                           boolean_array,probability,threshold,
                           no_event_parameters_object):
    # modify again, based on the integal noise
    x = split_fec.retract.Time
    x_sliced = x[slice_to_use]
    n_points = split_fec.tau_num_points
    force_sliced = split_fec.retract.Force[slice_to_use]
    boolean_ret = boolean_array.copy()
    probability_updated = probability.copy()
    interpolator = split_fec.retract_spline_interpolator(slice_to_use)
    interp_f = interpolator(x_sliced)
    probability_integral = _no_event.\
        _integral_probability(force_sliced,interp_f,n_points,
                              no_event_parameters_object)
    probability_updated[slice_to_use] *= probability_integral
    boolean_ret[slice_to_use] = (probability_updated[slice_to_use] < threshold)
    # zero out all the previous ones ie: no new points from this mask, just
    # better signal. XXX: not quite working right?
    where_not_already = np.where(np.logical_not(boolean_array))[0]    
    if (where_not_already.size > 0):
        probability_updated[where_not_already] = 1
        boolean_ret = (probability_updated < threshold)
    return slice_to_use,boolean_ret,probability_updated

def _condition_no_delta_significance(no_event_parameters_object,df_true,
                                     negative_only,interp_f,n_points):
    # XXX move to utility
    epsilon = no_event_parameters_object.epsilon
    sigma = no_event_parameters_object.sigma
    min_signal = (epsilon+sigma)
    epsilon_approach = no_event_parameters_object.delta_epsilon
    sigma_approach = no_event_parameters_object.delta_sigma
    if (negative_only):
        baseline = -min_signal
    else:
        # considering __all__ signal. XXX need absolute value df?
        baseline = min_signal
    if (negative_only):
        # XXX ?.... shouldnt this be minimum? (*dont* want positive)
        value_cond = (np.minimum(0,df_true) > baseline)
    else:
        # XXX should *not* need to have two separate methods. determine why
        # (probably adhesions)
        n_slice_region = df_true.size
        f0 = [interp_f[min(n_slice_region-1,i+n_points)] 
          for i in range(n_slice_region)]            
        interp_f_minus_baseline = interp_f - f0
        value_cond = (np.abs(interp_f_minus_baseline) < min_signal)
    return value_cond

def _condition_delta_at_zero(no_event_parameters_object,df_true,negative_only,
                             interp_f,pred_retract_surface_idx, slice_to_use):
    epsilon_approach = no_event_parameters_object.delta_epsilon
    sigma_approach = no_event_parameters_object.delta_sigma
    pred_retract_surface_idx_in_slice = pred_retract_surface_idx-\
                                        slice_to_use.start
    zero_force = interp_f[pred_retract_surface_idx_in_slice]
    diff = interp_f - np.maximum(0,df_true)
    zero_condition_baseline = zero_force+sigma_approach+epsilon_approach
    """
    plt.subplot(2,1,1)
    plt.plot(diff)
    plt.axhline(zero_force)
    plt.subplot(2,1,2)    
    plt.plot(diff)
    plt.axhline(zero_force)
    plt.show()
    """
    return (diff <= zero_condition_baseline) 



def delta_mask_function(split_fec,slice_to_use,
                        boolean_array,probability,threshold,
                        no_event_parameters_object,negative_only=True):
    x = split_fec.retract.Time
    force = split_fec.retract.Force
    x_sliced = x[slice_to_use]
    force_sliced = force[slice_to_use]
    n_points = split_fec.tau_num_points
    min_points_between = _min_points_between(n_points)
    boolean_ret = boolean_array.copy()
    probability_updated = probability.copy()
    # get the retract df spectrum
    interpolator = split_fec.retract_spline_interpolator(slice_to_use)
    interp_f = interpolator(x_sliced)
    df_true = _no_event._delta(x_sliced,interp_f,min_points_between)
    # get the baseline results
    ratio_probability = _no_event.\
        _delta_probability(df=df_true,
                           no_event_parameters=no_event_parameters_object,
                           negative_only=negative_only)
    probability_updated[slice_to_use] *= ratio_probability
    tol = 1e-9
    no_event_cond = (1-ratio_probability<tol)
    # find where the derivative is definitely not an event
    value_cond = \
        _condition_no_delta_significance(no_event_parameters_object,df_true,
                                         negative_only,interp_f,
                                         n_points)
    # find where we are consistent with zero
    pred_retract_surface_idx = split_fec.get_predicted_retract_surface_index()
    consistent_with_zero_cond = \
    _condition_delta_at_zero(no_event_parameters_object,df_true,negative_only,
                             interp_f,pred_retract_surface_idx,slice_to_use)
    gt_condition = np.ones(boolean_ret.size)
    gt_condition[slice_to_use] = ((value_cond) | (no_event_cond) | 
                                  (consistent_with_zero_cond))
    get_best_slice_func = lambda slice_list: \
        get_slice_by_max_value(interp_f,slice_to_use.start,slice_list)
    # update the boolean array before we slice
    boolean_ret = probability_updated < threshold
    boolean_ret,probability_updated = \
            safe_reslice(original_boolean=boolean_ret,
                         original_probability=probability_updated,
                         condition=gt_condition,
                         min_points_between=min_points_between,
                         get_best_slice_func=get_best_slice_func)
    boolean_ret = probability_updated < threshold
    """
    xlim = plt.xlim(min(x),max(x))
    plt.subplot(3,1,1)
    valid_idx = np.where(np.logical_not(gt_condition))
    invalid_idx = np.where(gt_condition)
    plt.plot(x[invalid_idx],force[invalid_idx],color='k',alpha=0.3)
    plt.plot(x[valid_idx],force[valid_idx],color='g')    
    plt.plot(x_sliced,interp_f,color='b')
    plt.xlim(xlim)
    plt.subplot(3,1,2)
    plt.plot(x_sliced,value_cond)
    plt.plot(x_sliced,no_event_cond+1.1)
    plt.xlim(xlim)
    plt.subplot(3,1,3)
    plt.semilogy(x,probability_updated,linestyle='--')
    plt.semilogy(x,probability)
    plt.xlim(xlim)
    plt.axhline(threshold)
    plt.xlim(xlim)
    plt.show()
    """
    return slice_to_use,boolean_ret,probability_updated

def adhesion_mask_function_for_split_fec(split_fec,slice_to_use,boolean_array,
                                         probability,threshold,
                                         no_event_parameters_object):
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
    n_points = split_fec.tau_num_points
    probability_functions_and_kw = \
        [ [derivative_mask_function,dict(negative_only=False)],
          [integral_mask_function,dict()],
          [delta_mask_function,dict(negative_only=False)]]
    probability_updated = probability.copy()      
    boolean_ret = boolean_array.copy()
    slice_updated = slice_to_use
    surface_index = split_fec.get_predicted_retract_surface_index()   
    for f,kw_tmp in probability_functions_and_kw:
        kw = dict(split_fec=split_fec,
                  slice_to_use=slice_updated,
                  threshold=threshold,
                  boolean_array=boolean_ret,
                  no_event_parameters_object=no_event_parameters_object,
                  probability=probability_updated,**kw_tmp)
        slice_update,boolean_ret,probability_tmp = f(**kw)
        probability_updated = probability_tmp*probability_updated 
    # determine where the surface is 
    non_events = probability_updated > threshold
    min_points_between = _min_points_between(n_points)
    min_idx = surface_index + min_points_between
    n = split_fec.retract.Force.size
    # remove all things before the predicted surface, and at the boundary
    boolean_ret[:min_idx] = 0
    boolean_ret[-min_points_between:] = 0
    probability_updated[:min_idx] = 1
    probability_updated[-min_points_between:] = 1
    slice_updated = slice(min_idx,n-min_points_between,1)
    # note: need to offset events
    no_event_mask = np.where(non_events)[0] 
    # XXX finish current event, keep consuming events until startd/end
    # are beyond threshold
    event_mask = np.where(~non_events)[0]
    if (event_mask.size ==0 or no_event_mask.size == 0):
        return slice_updated,boolean_ret,probability_updated
    # POST: we have at least one event and one non-event
    # (could be some adhesion!)
    # determine events that contain the surface index
    event_boundaries = _event_slices_from_mask(event_mask,min_points_between)
    # get a list of the events with a starting point below the surface
    events_containing_surface = [e for e in event_boundaries
                                 if (e.start <= surface_index)]     
    n_events_surface = len(events_containing_surface)
    if (n_events_surface == 0):
        return slice_updated,boolean_ret,probability_updated
    last_event_containing_surface_end = \
        events_containing_surface[-1].stop + min_points_between
    min_idx = max(min_idx,last_event_containing_surface_end)
    n = probability.size-1
    min_idx = min(n-(n_points+1),min_idx)       
    # update the boolean array and the probably to just reflect the slice
    # ie: ignore the non-unfolding probabilities above
    boolean_ret[:min_idx] = 0
    slice_updated = slice(min_idx,slice_updated.stop,1)
    """
    # set the interpolator for the non-adhesion region; need to re-calculate
    # since adhesion (probably) really screws everything up
    x = split_fec.retract.Time[slice_updated]
    y = split_fec.retract.Force[slice_updated]
    interp = split_fec.retract_spline_interpolator(slice_to_fit=slice_updated)
    split_fec.set_retract_knots(interp)
    # enable calculating the delta
    #no_event_parameters_object._set_valid_delta(True)
    # get the probability of only the negative regions
    probability_in_slice,_ = _no_event.\
        _no_event_probability(x,interp,y,n_points,no_event_parameters_object,
                              negative_only=True)
    """
    probability_updated = probability.copy()
    probability_updated[:min_idx] = 1
    return slice_updated,boolean_ret,probability_updated

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
    
def event_by_loading_rate(x,y,slice_event):
    """
    see _loading_rate_helper 

    Args:
        see _loading_rate_helper
    Returns:
        predicted index (absolute) in x,y where we think the event is happening
    """    
    fit_x,fit_y,pred,idx_above_predicted,local_max_idx = \
            _loading_rate_helper(x,y,slice_event)
    x_event = x[slice_event]
    y_event = y[slice_event]
    # POST: have a proper max, return the last time we are above
    # the linear prediction
    if (len(idx_above_predicted) == 0):
        return local_max_idx
    return idx_above_predicted[-1]

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
    interp_retract = split_fec.retract_spline_interpolator()
    # set the knots based on the initial interpolator, so that
    # any time we make a new splining object, we use the same knots
    split_fec.set_retract_knots(interp_retract)
    # set the epsilon and tau by the approach
    min_points_between = _min_points_between(n_points)    
    stdevs,epsilon,sigma,slice_fit_approach,spline_fit_approach =\
        split_fec._approach_metrics()
    split_fec.set_espilon_and_sigma(epsilon,sigma)
    split_fec.set_approach_metrics(slice_fit_approach,spline_fit_approach)
    """
    get the interpolator delta in the slice
    """
    interpolator_approach_x = split_fec.approach.Time[slice_fit_approach]
    interpolator_approach_f = spline_fit_approach(interpolator_approach_x)
    df_approach = Analysis.local_centered_diff(interpolator_approach_f,
                                               n=min_points_between)
    delta_epsilon,delta_sigma = np.median(df_approach),np.std(df_approach)
    """
    get the interpolate derivative in the slice
    """
    approach_interp_deriv = \
            spline_fit_approach.derivative()(interpolator_approach_x)
    derivative_epsilon = np.median(approach_interp_deriv)
    derivative_sigma = np.std(approach_interp_deriv)
    # get the remainder of the approach metrics needed
    # note: to start, we do *not* use delta; this is calculated
    # after the adhesion
    approach_dict = dict(integral_sigma   = 2*sigma*min_points_between,
                         integral_epsilon = epsilon,
                         delta_epsilon = delta_epsilon,
                         delta_sigma   = delta_sigma,
                         derivative_epsilon = derivative_epsilon,
                         derivative_sigma   = derivative_sigma,
                         valid_delta = False,
                         valid_integral = False,
                         valid_derivative = False,
                         **kwargs)
    # call the predict function
    final_kwargs = dict(epsilon=epsilon,sigma=sigma,**approach_dict)
    to_ret = _predict(x=time,
                      y=force,
                      n_points=n_points,
                      interp=interp_retract,
                      threshold=threshold,
                      local_event_idx_function=event_by_loading_rate,
                      **final_kwargs)
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
    f_refs = [adhesion_mask_function_for_split_fec,delta_mask_function]
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

