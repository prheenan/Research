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

def local_stdev(f,n):
    """
    Gets the local standard deviaiton (+/- n, except at boundaries 
    where it is just in the direction with data

    Args:
        f: what we want the stdev of
        n: window size
    Returns:
        array, same size as f, with the dat we want
    """
    max_n = f.size
    return [np.std(f[max(0,i-n):min(max_n,i+n)]) for i in range(max_n)]

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
    time,separation,force = retract.Time,retract.Separation,retract.Force
    n_points = ex.tau_num_points
    interp = ex.retract_spline_interpolator(deg=2)
    interpolated_force = interp(time)
    event_slices = ex.get_retract_event_idx()
    # get a model for the local standard deviaiton using the autocorrelation
    # time from the event
    diff = force-interpolated_force
    stdevs = local_stdev(diff,n_points)
    # plot everything
    style_events = dict(color='r',label="True events")
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    plt.plot(time,force,color='k',alpha=0.3)
    plt.plot(time,interpolated_force,color='b',linewidth=2)
    Plotting.highlight_events(event_slices,time,force,
                              **style_events)
    PlotUtilities.lazyLabel("","Force (au)","")
    plt.subplot(2,1,2)
    plt.plot(time,stdevs)
    # plot the autocorrelation time along the plot
    min_x_auto = min(time) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+ex.tau]
    plt.plot(auto_correlation_x, [np.max(stdevs),np.max(stdevs)],
             linewidth=5,color='g',label="autocorrelation time")
    Plotting.highlight_events(event_slices,time,stdevs,linewidth=5,
                              **style_events)
    PlotUtilities.lazyLabel("Time (au)","Event Stdev","")
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
