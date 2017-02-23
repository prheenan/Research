# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from Research.Personal.EventDetection.OtherMethods import method_helper
from Research.Personal.EventDetection.Util import Analysis,Plotting,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector

from GeneralUtil.python import PlotUtilities

def plot_prediction_info(ex,info,xlabel="Time",
                         x_func=lambda ex: ex.Time -ex.Time[0]):   
    """
    plots prediction information for the no-hypothesis predictor

    Args:
        ex: instance of (zeroed, etc) TimeSepForce instance
    
    Returns:
        This is a description of what is returned.
    """
    thresh = info.threshold
    retract = ex.retract
    event_slices = ex.get_retract_event_idx()
    time,separation,force = retract.Time,retract.Separation,retract.Force
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
    n_plots = 3
    x = x_func(retract)
    min_x,max_x = min(x),max(x)
    ylabel = "Force [pN]"
    x_range = max_x - min_x
    fudge = x_range * 0.05
    x_limits = [min_x - fudge,max_x + fudge]
    force_plot = force * 1e12
    plt.subplot(n_plots,1,1)
    plt.plot(x,force_plot,color='k',alpha=0.3)
    plt.plot(x,interpolated_force,color='b',linewidth=2)
    Plotting.highlight_events(event_slices,x,force_plot,
                              **style_events)
    PlotUtilities.lazyLabel("",ylabel,"")
    plt.xlim(x_limits)
    plt.subplot(n_plots,1,2)
    # plot the autocorrelation time along the plot
    min_x_auto = min(x) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+tau]
    plt.semilogy(x[info.slice_fit],info.cdf)
    Plotting.highlight_events(event_slices,x,info.cdf,linewidth=5,
                              **style_events)
    plt.axhline(thresh,label="threshold",linestyle='--',color='r')
    PlotUtilities.lazyLabel("","No-Event CDF ","")
    plt.xlim(x_limits)
    plt.subplot(n_plots,1,3)
    # XXX check mask has at least one...
    plt.plot(x,force_plot,'b-',color='k',alpha=0.3)
    for fwd,rev,event in zip(event_idx_end,event_idx_start,event_idx):
        plt.axvline(x[fwd],linestyle='--',color='r')
        plt.axvline(x[rev],color='g')
        plt.axvline(x[event],linewidth=3)
    plt.xlim(x_limits)
    PlotUtilities.lazyLabel(xlabel,ylabel,"")

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    ex = method_helper.get_example()
    thresh = 5e-3
    info = Detector._predict_helper(ex,threshold=thresh)
    # XXX fix threshhold
    fig = PlotUtilities.figure(figsize=(8,12))    
    plot_prediction_info(ex,info)
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
