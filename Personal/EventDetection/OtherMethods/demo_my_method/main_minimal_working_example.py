# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from Research.Personal.EventDetection.OtherMethods import method_helper
from Research.Personal.EventDetection.Util import Analysis,Plotting,Scoring

from GeneralUtil.python import PlotUtilities
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
    def __init__(self,event_idx,start_idx,end_idx,local_stdev,interp,mask,
                 cdf,slice_fit):
        """
        record the data from _predict_helper

        Args:
            event_idx : the event centers we found
            start_idx : the start of the events
            end_idx  : the end of the events
            local_stdev : the local standad deviation at each point in slice_fit
            interp_mask : points considered valid in slice_fit
            cdf : cummulative density function for the probability of a point
            in slice_fit, given the model of the data we have
        
            slice_fit : slice within the original retract where we tried to 
            find events. We have to remove the first few and last few points
        Returns:
            prediction_info object
        """
        self.event_idx = event_idx
        self.start = start_idx
        self.end = end_idx
        self.local_stdev = local_stdev
        self.interp = interp
        self.mask = mask
        self.cdf = cdf
        self.slice_fit = slice_fit

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
    min_points_between = int(np.ceil(n_points/2))
    # N degree b-spline has continuous (N-1) derivative
    interp = split_fec.retract_spline_interpolator(deg=3)
    interp_first_deriv = interp.derivative(1)(time)
    # get the interpolated derivative
    interpolated_force = interp(time)
    # get a model for the local standard deviaiton using the autocorrelation
    # time from the event
    diff = force-interpolated_force
    stdevs = local_stdev(diff,n_points)
    # get the cwt of the wavelet; see pp219 of Mallat, Wavelet Tour (XXX TODO)
    global_stdev = np.std(diff)
    median_local_stdev = np.median(stdevs)
    slice_fit = slice(min_points_between,-min_points_between,1)
    stdev_masked = stdevs[slice_fit]
    q25,q75 = np.percentile(stdev_masked,[25,75])
    iqr = q75-q25
    cdfs = 1-stats.norm.cdf(stdev_masked,loc=median_local_stdev,scale=iqr)
    mask = np.where(cdfs <=  threshold)[0]
    # add back in the offset
    mask += min_points_between
    last_point = mask[-1]
    idx_step_changes = np.where(np.diff(mask) > min_points_between)[0]
    step_end_idx = mask[idx_step_changes]
    step_start_idx = mask[(idx_step_changes + 1)]
    event_idx_end = list(step_end_idx) + [mask[-1]] 
    event_idx_start = [mask[0]] + list(step_start_idx)
    event_slices = [slice(start,end,1) 
                    for start,end in zip(event_idx_start,event_idx_end)]
    # determine where the first derivative is minimal (most negative, XXX check)
    # in each slice; that is the strongest indicator that an event is taking 
    # place
    min_deriv_idx = [e.start + np.argmin(interp_first_deriv[e])
                     for e in event_slices]
    event_idx = min_deriv_idx
    to_ret = prediction_info(event_idx = event_idx,
                             start_idx = event_idx_start,
                             end_idx   = event_idx_end,
                             local_stdev = stdevs,
                             interp = interp,
                             mask = mask,
                             cdf=cdfs,
                             slice_fit=slice_fit)
    return to_ret

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    ex = method_helper.get_example()
    retract = ex.retract
    event_slices = ex.get_retract_event_idx()
    time,separation,force = retract.Time,retract.Separation,retract.Force
    # XXX fix threshhold
    thresh = 1e-3
    info = _predict_helper(ex,threshold=thresh)
    event_idx_end,event_idx_start,event_idx = info.end,\
                                              info.start,\
                                              info.event_idx
    mask = info.mask
    interp_first_deriv = info.interp.derivative(1)(time)
    # get the interpolated derivative
    interpolated_force = info.interp(time)
    tau = ex.tau
    stdevs = info.local_stdev
    # plot everything
    style_events = dict(color='r',label="True events")
    fig = PlotUtilities.figure()
    n_plots = 3
    x = separation
    min_x,max_x = min(x),max(x)
    x_range = max_x - min_x
    fudge = x_range * 0.05
    x_limits = [min_x - fudge,max_x + fudge]
    x_label = "Time"
    plt.subplot(n_plots,1,1)
    plt.plot(x,force,color='k',alpha=0.3)
    plt.plot(x,interpolated_force,color='b',linewidth=2)
    Plotting.highlight_events(event_slices,time,force,
                              **style_events)
    PlotUtilities.lazyLabel("","Force (au)","")
    plt.xlim(x_limits)
    plt.subplot(n_plots,1,2)
    # plot the autocorrelation time along the plot
    min_x_auto = min(x) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+ex.tau]
    plt.semilogy(x[info.slice_fit],info.cdf)
    Plotting.highlight_events(event_slices,x,info.cdf,linewidth=5,
                              **style_events)
    plt.axhline(thresh,label="threshold",linestyle='--',color='r')
    PlotUtilities.lazyLabel("","CDF ","")
    plt.xlim(x_limits)
    plt.subplot(n_plots,1,3)
    # XXX check mask has at least one...
    plt.plot(x,force,'b-',color='k',alpha=0.3)
    for fwd,rev,event in zip(event_idx_end,event_idx_start,event_idx):
        plt.axvline(x[fwd],linestyle='--',color='r')
        plt.axvline(x[rev],color='g')
        plt.axvline(x[event],linewidth=3)
    plt.xlim(x_limits)
    PlotUtilities.lazyLabel(x_label,"Force","")
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
