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
from scipy import signal

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
    def __init__(self,event_idx,start_idx,end_idx,local_stdev,interp,mask):
        self.event_idx = event_idx
        self.start = start_idx
        self.end = end_idx
        self.local_stdev = local_stdev
        self.interp = interp
        self.mask = mask

def _predict_helper(split_fec,threshold):
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
    mask = np.where(stdevs[min_points_between:] >  threshold)[0]
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
                             mask = mask)
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
    thresh = 0.6e-11
    info = _predict_helper(ex,threshold=thresh)
    event_idx_end,event_idx_start,event_idx = info.end,\
                                              info.start,\
                                              info.event_idx
    print(event_idx)
    mask = info.mask
    interp_first_deriv = info.interp.derivative(1)(time)
    # get the interpolated derivative
    interpolated_force = info.interp(time)
    tau = ex.tau
    stdevs = info.local_stdev
    # plot everything
    style_events = dict(color='r',label="True events")
    time_limits = [min(time),max(time)]
    fig = PlotUtilities.figure()
    n_plots = 4
    plt.subplot(n_plots,1,1)
    plt.plot(time,force,color='k',alpha=0.3)
    plt.plot(time,interpolated_force,color='b',linewidth=2)
    Plotting.highlight_events(event_slices,time,force,
                              **style_events)
    PlotUtilities.lazyLabel("","Force (au)","")
    plt.xlim(time_limits)
    plt.subplot(n_plots,1,2)
    plt.plot(time,stdevs)
    # plot the autocorrelation time along the plot
    min_x_auto = min(time) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+ex.tau]
    plt.plot(auto_correlation_x, [np.max(stdevs),np.max(stdevs)],
             linewidth=5,color='g',label="autocorrelation time")
    Plotting.highlight_events(event_slices,time,stdevs,linewidth=5,
                              **style_events)
    PlotUtilities.lazyLabel("","Local Stdev","")
    plt.axhline(thresh,label="threshold",linestyle='--',color='r')
    plt.xlim(time_limits)
    plt.subplot(n_plots,1,3)
    # XXX check mask has at least one...
    plt.plot(time,stdevs,'b.')
    for fwd,rev,event in zip(event_idx_end,event_idx_start,event_idx):
        plt.axvline(time[fwd],linestyle='--')
        plt.axvline(time[rev])
        plt.axvline(time[event],linewidth=3)
    plt.xlim(time_limits)
    PlotUtilities.lazyLabel("","","")
    plt.subplot(n_plots,1,4)
    plt.plot(time[mask],interp_first_deriv[mask],'b.')
    plt.xlim(time_limits)
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
