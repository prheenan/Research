# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import signal,stats

from Research.Personal.EventDetection.Util import Analysis,Plotting
from GeneralUtil.python import PlotUtilities,GenUtilities
from itertools import chain
from scipy.signal import medfilt

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

def safe_cheby_probability(y,loc,scale):
    to_ret = np.ones(y)
    k = (possible_deriv-loc)/scale
    possible_idx = np.where(k > 1)
    to_ret[possible_idx] = 1/k[possible_idx]**2
    return to_ret

def _spline_derivative_probability_generic(x,interpolator,scale=None,loc=None):
    """
    see  spline_derivative_probability, except a genertic method
    
    Args:
        x: x values
        interpolator: to interpolate along
    Returns:
        see spline_derivative_probability
    """ 
    derivative_force = interpolator.derivative()(x)
    if (loc is None):
        loc = np.median(derivative_force)
    if (scale is None):
        std_iqr = np.std(derivative_force)
        scale = std_iqr
    probability  = _no_event_chebyshev(derivative_force,loc,scale)
    return probability 

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
    median = np.median(interp_f)
    bool_interp = ( (interp_f - stdev < median) |
                    (stdev - epsilon < sigma) )
    where_not_bool_in_slice = np.where(~bool_interp)[0]
    no_event_possible = np.ones(boolean_array.size)
    """
    # XXX debugging
    plt.subplot(2,1,1)
    plt.plot(f,color='k',alpha=0.3)
    plt.plot(x,interp_f)
    plt.plot(x[where_not_bool_in_slice],interp_f[where_not_bool_in_slice],
             color='r')
    plt.subplot(2,1,2)
    plt.plot(stdev-epsilon)
    plt.axhline(thresh,linestyle='--')
    plt.axhline(-thresh,linestyle='--')
    plt.show()
    """
    no_event_possible[slice_to_use] = bool_interp
    get_best_slice_func = lambda slice_list: \
        get_slice_by_max_value(interp_f,slice_to_use.start,slice_list)
    ret = safe_reslice(boolean_ret,probability_updated,
                       condition=no_event_possible,
                       min_points_between=min_points_between,
                       get_best_slice_func=get_best_slice_func)
    boolean_updated,probability_updated = ret
    """
    # XXX debugging...
    Plotting.debug_plot_force_value(x,f,interp_f,probability,
                                    probability_updated,
                                    slice_to_use,bool_interp)
    plt.show()
    """
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
    """
    # XXX debugging
    ylim = [min(original_probability)/2,2]
    plt.subplot(3,1,1)
    plt.plot(condition,label="condition")
    idx = np.arange(condition.size)
    colors = ['r','g']
    for i,e in enumerate(new_events):
        plt.plot(idx[e],condition[e],color=colors[i % 2],linewidth=3)
    plt.ylim([-0.2,1.1])
    PlotUtilities.lazyLabel("","","")
    plt.subplot(3,1,2)
    plt.plot(new_boolean+1.1,label="new")
    plt.plot(new_probability)
    plt.ylim(ylim)
    PlotUtilities.lazyLabel("index","","")
    plt.yscale("log")
    plt.subplot(3,1,3)
    plt.plot(original_boolean+1.1,label="original")
    plt.plot(original_probability)
    plt.yscale("log")
    plt.ylim(ylim)
    PlotUtilities.lazyLabel("index","","")
    plt.show()
    """
    return new_boolean,new_probability

def derivative_mask_function(split_fec,slice_to_use,
                             boolean_array,probability,threshold,
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
    if ((slice_to_use.stop - offset) < min_points_between):
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
    offset = slice_to_use.start    
    force_sliced = force[slice_to_use]
    interp_sliced = interp(x_sliced)
    interp_slice_deriv = interp(x_sliced,1)
    diff_sliced = interp_sliced - force_sliced
    local_stdev = Analysis.local_stdev(diff_sliced,n_points)
    epsilon,sigma = split_fec.get_epsilon_and_sigma()
    ratio = (interp_slice_deriv*split_fec.tau)/local_stdev
    ratio_min_threshold = -1
    # XXX debuugging...
    idx_offset_approach = split_fec.get_predicted_approach_surface_index()
    n_approach = split_fec.approach.Force.size
    offset_approach = max([offset,n_approach - idx_offset_approach,
                           n_approach-absolute_max_index])
    slice_approach = slice(offset_approach,-offset_approach,1)
    approach_interp = split_fec.\
        approach_spline_interpolator(slice_to_fit=slice_approach)
    approach = split_fec.approach
    approach_force = approach.Force[slice_approach]
    approach_time = approach.Time[slice_approach]
    approach_interp_sliced = approach_interp(approach_time)
    approach_interp_deriv = approach_interp(approach_time,1)
    min_deriv = np.min(approach_interp_deriv)
    med_deriv_appr = np.median(approach_interp_deriv)
    std_deriv_appr = np.std(approach_interp_deriv)
    kwargs_approach_deriv = dict(loc=med_deriv_appr,
                                 scale=std_deriv_appr)
    # modulate the probabilities by the approach
    """
    #XXX debugging
    Plotting.debug_plot_derivs(approach_time,approach_force,
                      approach_interp_sliced,x_sliced,
                      force_sliced,interp_sliced,
                      approach_interp_deriv,interp_slice_deriv,
                      min_deriv)
    plt.show()
    """
    prob_mod = _spline_derivative_probability_generic(x_sliced,interp,
                                                      **kwargs_approach_deriv)
    probability_updated[slice_to_use] *= prob_mod
    # modify again, based on the integal noise
    interpolator = split_fec.retract_spline_interpolator(slice_to_use)
    interp_f = interpolator(x_sliced)
    diff = force_sliced-interp_f
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
    """
    #XXX debugging
    xlim = [min(time),max(time)]
    plt.subplot(4,1,1)
    plt.plot(x_sliced,force_sliced,alpha=0.3,color='k')
    plt.plot(x_sliced,interp_f)
    plt.xlim(xlim)
    plt.subplot(4,1,2)
    plt.semilogy(x_sliced,prob_mod,label='new',linestyle='--')
    plt.semilogy(x_sliced,probability[slice_to_use]*prob_mod,color='r',
                label='old')
    plt.xlim(xlim)  
    plt.legend()    
    plt.subplot(4,1,3)
    plt.plot(x_sliced,local_integral)
    plt.axhline(sigma+epsilon,color='r',linestyle='--')
    plt.axhline(thresh_integral)
    plt.xlim(xlim)  
  
    plt.subplot(4,1,4)
    plt.semilogy(time,probability_updated,color='g',label='updated',
                 linestyle='--')
    plt.semilogy(time,probability,color='r',label='old')
    plt.xlim(xlim)
    plt.legend()
    plt.show()
    # XXX debugging
    xlim = [min(time),max(time)]
    plt.subplot(2,1,1)
    plt.plot(x_sliced,interp(x_sliced))
    plt.xlim(xlim)
    plt.subplot(2,1,2)
    plt.plot(time,probability_updated)
    plt.xlim(xlim)
    plt.yscale("log")
    plt.show()
    """
    boolean_ret[slice_to_use] = (probability_updated[slice_to_use] < threshold)
    # find where the derivative is definitely not an event
    gt_condition = np.ones(boolean_ret.size)
    """
    # XXX debugging filtering
    median_filter_points = min_points_between
    if (median_filter_points % 2 == 0):
        median_filter_points += 1
    median_filtered = medfilt(interp_sliced,median_filter_points)
    plt.plot(interp_sliced-median_filtered)
    plt.axhline(sigma)
    plt.show()
    """
    gt_condition[slice_to_use] = ( (interp_slice_deriv > 0) |
                                   (interp_sliced - stdev < median) |
                                   (ratio > ratio_min_threshold))
    get_best_slice_func = lambda slice_list: \
        get_slice_by_max_value(interp_sliced,slice_to_use.start,slice_list)
    boolean_ret,probability_updated = \
            safe_reslice(original_boolean=boolean_ret,
                         original_probability=probability_updated,
                         condition=gt_condition,
                         min_points_between=min_points_between,
                         get_best_slice_func=get_best_slice_func)
    """
    #XXX debugging
    Plotting.debug_plot_derivative_ratio(time,slice_to_use,
                                         ratio,interp_sliced,force_sliced,
                                         interp_slice_deriv,
                                         boolean_ret,probability_updated,
                                         absolute_min_idx,ratio_min_threshold)
    plt.show()
    #XXX debugging...
    xlim = [min(time),max(time)]
    plt.subplot(2,1,1)
    plt.plot(x_sliced,force_sliced*1e12)
    plt.xlim(xlim)
    PlotUtilities.lazyLabel("","Force","")
    plt.subplot(2,1,2)
    plt.semilogy(time,probability_updated)
    plt.plot(time,boolean_array + 1.1,label="original")
    plt.plot(time,boolean_ret + 1.1,label="updated")
    plt.xlim(xlim)
    PlotUtilities.lazyLabel("Time","masks","")
    plt.show()
    """
    return slice_to_use,boolean_ret,probability_updated

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
    probability_updated= probability.copy()
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
        epsilon,sigma = split_fec.get_epsilon_and_sigma()
        kwargs = dict(epsilon=epsilon,sigma=sigma)
        probability_updated_slice, _ = \
            _no_event_probability(time,interp,force,n_points,
                                  slice_fit=slice_update,**kwargs)
        probability_updated[slice_update] = probability_updated_slice
    else:
        slice_update = slice(0,None,1)        
    return slice_update,boolean_array,probability_updated

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
    retract = split_fec.retract   
    # determine the boundaries of the 'no events'
    min_points_between = _min_points_between(n_points)
    min_idx = surface_index + min_points_between
    # remove all things before the predicted surface, and at the boundary
    to_ret[:min_idx] = 0    
    to_ret[-min_points_between:] = 0    
    # determine where the derivative after the predicted surface is first zero
    time = retract.Time
    # we *fit* the spline to the entire data set
    slice_to_fit = slice(0,None,1)
    # we only *consider* points that are above some marker (e.g. surface_index)
    slice_idx = np.arange(0,time.size)
    time_to_fit = time[slice_to_fit]
    force_to_fit = retract.Force[slice_to_fit]
    interp = split_fec.retract_spline_interpolator(slice_to_fit=slice_to_fit)
    force_fit = interp(time_to_fit)
    deriv = interp.derivative()(time_to_fit)
    derivative_gt_zero = deriv > 0
    derivative_le_zero = deriv <= 0
    # must be above the surface *and* have a derivative <= 0 
    where_deriv_le_zero = \
        np.where(derivative_le_zero & (slice_idx > surface_index) )[0]
    if (where_deriv_le_zero.size == 0):
        return to_ret
    # POST: have at least one point we can test, figure out where we get back
    # to the force baseline
    where_derivative_le_zero_absolute = where_deriv_le_zero
    min_idx = where_derivative_le_zero_absolute[0]
    to_ret[:min_idx] = 0
    # POST: we found a peak (or a flat point) followed by some kind of increase
    # this means we should have passed (at least one) adhesion peak.
    # however, if there is an event happening here, we need to continue
    # consuming them (noting that we determine the probability distribution
    # again, since that initial adhesion prolly messed us up 
    force = retract.Force
    epsilon,sigma = split_fec.get_epsilon_and_sigma()
    kwargs_no_event = dict(epsilon=epsilon,sigma=sigma)
    probability_distribution = _no_event_probability(x=time,
                                                     interp=interp,
                                                     y=force,
                                                     n_points=n_points,
                                                     **kwargs_no_event)
    non_events = probability_distribution > threshold    
    # note: need to offset events
    no_event_mask = np.where(non_events)[0] + min_idx
    # XXX finish current event, keep consuming events until startd/end
    # are beyond threshold
    event_mask = np.where(~non_events)[0] + min_idx
    if (event_mask.size ==0 or no_event_mask.size == 0):
        return to_ret
    # POST: we have at least one event and one non-event
    # (could be some adhesion!)
    # determine events that contain the surface index
    event_boundaries = _event_slices_from_mask(event_mask,min_points_between)
    # get a list of the events with a starting point below the surface
    events_containing_surface = [e for e in event_boundaries
                                 if (e.start <= min_idx)]
    n_events_surface = len(events_containing_surface)  
    if (n_events_surface == 0):
        return to_ret
    # POST: at least one event contains the surface. Update the minimum index
    # to go to the end of the (last) event below or at the surface, unless
    # the end is below the surface, then just keep the minimum at the surface
    last_event_containing_surface_end = \
        events_containing_surface[-1].stop + min_points_between
    min_idx = max(min_idx,last_event_containing_surface_end)    
    to_ret[:min_idx] = 0              
    # finally, make sure we are at least at the N-th percentile
    N = 75
    pct_threshold = np.percentile(force_to_fit[min_idx:],N)
    where_below = np.where((force_to_fit < pct_threshold) & \
                           (slice_idx > min_idx))[0]
    if (where_below.size > 0):
        min_idx = where_below[0]
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

def _no_event_chebyshev(g,epsilon,sigma):
    denom = (g-epsilon)
    k_chebyshev = denom/sigma
    # determine where the chebyshev is 'safe', otherwise we are at or above
    # the mean estimate and hence not a useful metric
    cheby_idx = np.where(k_chebyshev >= 1)
    chebyshev = np.ones(k_chebyshev.size)
    # actually calculate the upper bound for the probability
    chebyshev[cheby_idx] = (1/k_chebyshev[cheby_idx])**2
    return chebyshev

def _no_event_probability(x,interp,y,n_points,epsilon=None,sigma=None,
                          slice_fit=None):
    """
    returns the no-event probability at each point in y

    Args:
        x: the x values that interp takes, see _event_probabilities
        y: the y values we are searching for an event, see _event_probabilities
        interp: see _event_probabilities
        n_points: number of points to use in estimating r(q)=g-g* by the 
        local standard deviaiton of y-interp(x)
        epsilon: the mean error parameter
        sigma: the stdev of the error parameter
        slice_fit: an optional slice to use to compute the probabilities
    Returns:
        tuple of <probability, local stdevs>
    """
    if (slice_fit is None):
        slice_fit = slice(0,None,1)
    n_original = x.size
    x = x[slice_fit]
    y = y[slice_fit]
    # get the interpolated function
    interpolated_y = interp(x)
    stdev_masked,epsilon_def,sigma_def = Analysis.\
        stdevs_epsilon_sigma(y,interpolated_y,n_points)
    if epsilon is None:
        epsilon = epsilon_def
    if sigma is None:
        sigma = sigma_def
    # note: chebyshev is like
    # P(|X - mu| >=  k * sigma) <= 1/k^2
    # we write k = (s(q) - epsilon)/scale
    chebyshev = _no_event_chebyshev(stdev_masked,epsilon,sigma)
    # for the edge cases, assume the probability is one                         
    probability_distribution = np.ones(y.size)
    # get the probability for all the non edge cases
    probability_distribution = chebyshev
    return probability_distribution,stdev_masked
        
def _event_probabilities(x,y,interp,n_points,threshold,**kwargs):
    """
    determines the mask (and associated event detection information)
    
    Args:
        x,y: independent and dependent variable (ie: 'q' and 'g'
        interp: the approximation to y vs x (ie: g*)
        n_points: number of points from autocorrelation function (ie: tau)

        threshold: maximum probability that a given datapoint fits the 
        model
    
        **kwargs; passed to no_event_probability
    Returns:
        tuple of :
            probability_distribution : no-event probability for each point in y
            slice_fit : the part of x and y that mask is valid for
            stdevs: the local, windowed standard deviation, s(q)
    """
    min_points_between = _min_points_between(n_points)
    slice_fit = slice(min_points_between,-min_points_between,1)
    probability_distribution = np.ones(x.size)
    probability_distribution_slice,stdevs = \
        _no_event_probability(x,interp,y,n_points=n_points,slice_fit=slice_fit,
                              **kwargs)
    probability_distribution[slice_fit] = probability_distribution_slice 
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


def _predict(x,y,n_points,interp,threshold,local_event_idx_function,
             remasking_functions=None,**kwargs):
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
    
        remasking_functions: funcitons that take in the slice, boolean array, 
        probability, and threshold and return update values for all except
        threshold
    
        kwargs: passed to _event_probabilities
    Returns:
        list of event slices
    """
    min_points_between = _min_points_between(n_points)
    probability_distribution,slice_fit,stdevs = \
        _event_probabilities(x,y,interp,n_points,threshold,**kwargs)
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
            # mask on probability distribution, to keep things consistent
            probabilities.append(probability_distribution)
            masks.append(np.where(bool_array)[0])
    # only keep points where we are farther than min_points between from the 
    # edges (ie: from index 0 and N-1)
    mask = np.where(bool_array)[0]
    n = mask.size
    if (mask.size > 0):
        event_slices = _event_slices_from_mask(mask,min_points_between)
    else:
        event_slices = []
    # XXX reject events with a very small time?
    event_duration = [ (e.stop-e.start) for e in event_slices]
    delta_split_rem = [ int(np.ceil(min_points_between-(delta))/2)
                        for delta in event_duration]
    # determine where the events are happening locally (guarentee at least
    # a search window of min_points)
    remainder_split = [max(0,d) for d in delta_split_rem ]
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
    interp = split_fec.retract_spline_interpolator()
    # set the knots based on the initial interpolator, so that
    # any time we make a new splining object, we use the same knots
    split_fec.set_retract_knots(interp)
    # set the epsilon and tau by the approach
    min_points_between = _min_points_between(n_points)    
    epsilon,sigma = split_fec.calculate_epsilon_and_sigma(n_points=n_points)
    split_fec.set_espilon_and_sigma(epsilon,sigma)
    final_kwargs = dict(epsilon=epsilon,sigma=sigma,**kwargs)
    to_ret = _predict(x=time,
                      y=force,
                      n_points=n_points,
                      interp=interp,
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

