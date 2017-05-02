# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting,Learning
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from mpl_toolkits.axes_grid.inset_locator import inset_axes,mark_inset

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    data_file = "../_Data/example_protein.pkl"
    data = CheckpointUtilities.lazy_load(data_file)
    split_fec,info_final = Detector._predict_full(data)
    # get the plotting versions of the time, etc
    time= split_fec.retract.Time
    time_plot = time - time[0]
    force_plot = split_fec.retract.Force * 1e12
    force_interp_plot = split_fec.retract_spline_interpolator()(time) * 1e12
    # plot everything
    raw_force_kwargs = dict(color='k',alpha=0.3)
    interp_force_kwargs = dict(color='b',linewidth=3)
    probabiity_kwargs = dict(color='r')
    n_cols = 2
    n_rows = 2
    ylim_force_pN = [-30,max(force_interp_plot)*1.1]
    ylim_prob = [min(info_final.cdf)/2,2]
    fig = PlotUtilities.figure(figsize=(8,8))
    # plot the 'raw' force
    ax1 = plt.subplot(n_rows,n_cols,1)
    plt.plot(time_plot,force_plot)
    plt.plot(time_plot,force_interp_plot)
    PlotUtilities.lazyLabel("","Force (pN)","")
    plt.ylim(ylim_force_pN)
    PlotUtilities.no_x_label()
    ax_zoom = plt.subplot(n_rows,n_cols,3)
    zoom_slice = info_final.event_slices[1]
    time_slice = time_plot[zoom_slice]
    plt.plot(time_slice,force_plot[zoom_slice])
    plt.plot(time_slice,force_interp_plot[zoom_slice])
    xlim = [time_slice[0],time_slice[-1]]
    plt.xlim(xlim)
    PlotUtilities.no_x_label()
    PlotUtilities.zoom_effect01(ax1, ax_zoom, *xlim)
    PlotUtilities.lazyLabel("Time (s)","Force (pN)","")    
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
