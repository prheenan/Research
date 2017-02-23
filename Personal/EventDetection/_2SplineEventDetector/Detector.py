# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import signal,stats


def local_stdev(f,n):
    """
    Gets the local standard deviaiton (+/- n, except at boundaries 
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

class prediction_info:
    def __init__(self,event_idx,event_slices,local_stdev,interp,mask,
                 cdf,slice_fit,threshold):
        """
        record the data from _predict_helper

        Args:
            event_idx : the event centers we found
            start_idx : the boundaries found for events
            local_stdev : the local standad deviation at each point in slice_fit
            interp_mask : points considered valid in slice_fit
            cdf : cummulative density function for the probability of a point
            in slice_fit, given the model of the data we have
        
            slice_fit : slice within the original retract where we tried to 
            find events. We have to remove the first few and last few points
            
            threshold: the threshhold for determining an event
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

def _event_mask(probability,threshold):
    """
    Given a probability distribution and a threshold, returns the indices
    where the probability is less  than the threshhold
    
    Args;
        probability: array of numbers between 0 and 1
        threshold: maximum value in mask
    Returns:
        indices where the probability is less than the threshold
    """
    return np.where(probability <= threshold)[0]
        
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
    min_points_between = int(np.ceil(n_points/2))    
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
    k_chebyshev = (stdev_masked-median_local_stdev)/scale
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
    mask_offset = mask[min_points_between:-min_points_between]
    idx_step_changes = np.where(np.diff(mask_offset) >= min_points_between)[0]
    step_end_idx = mask_offset[idx_step_changes]
    step_start_idx = mask_offset[(idx_step_changes + 1)]
    event_idx_end = list(step_end_idx) + [mask_offset[-1]] 
    event_idx_start = [mask_offset[0]] + list(step_start_idx)
    offset = min_points_between
    event_slices = [slice(start+offset,end+offset,1) 
                    for start,end in zip(event_idx_start,event_idx_end)]    
    return event_slices
    
def _predict_helper(split_fec,threshold):
    """
    uses spline interpolation and local stadard deviations to predict
    events.

    Args:
        split_fec: split_force_extension object, already initialized, and 
        zerod, with the autocorraltion time set. 

        threshhold: maximum probability that a given datapoint fits the 
        model
    Returns:
        prediction_info object
    """
    retract = split_fec.retract
    time,separation,force = retract.Time,retract.Separation,retract.Force
    n_points = split_fec.tau_num_points
    # N degree b-spline has continuous (N-1) derivative
    interp = split_fec.retract_spline_interpolator(deg=2)
    min_points_between = int(np.ceil(n_points/2))    
    probability_distribution,slice_fit,stdevs = \
        _event_probabilities(time,force,interp,n_points,threshold)
    mask = _event_mask(probability_distribution,threshold)
    if (mask.size > 0):
        event_slices = _event_slices_from_mask(mask,min_points_between)
    else:
        event_slices = []
    # determine where the first derivative is minimal (most negative, XXX check)
    # in each slice; that is the strongest indicator that an event is taking 
    # place
    max_force_idx = [e.start + np.argmax(force[e]) for e in event_slices]
    # XXX probably want to walk back up to the maximum force?
    event_idx = max_force_idx
    to_ret = prediction_info(event_idx = event_idx,
                             event_slices = event_slices,
                             local_stdev = stdevs,
                             interp = interp,
                             mask = mask,
                             cdf=probability_distribution,
                             slice_fit=slice_fit,
                             threshold=threshold)
    return to_ret
