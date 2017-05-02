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

from Research.Personal.EventDetection.Util import Offline,Plotting,Learning,\
    Analysis
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
    ylim_force_pN = [-40,max(force_interp_plot)*1.2]
    ylim_prob = [min(info_final.cdf)/2,2]
    fig = PlotUtilities.figure(figsize=(8,8))
    # # plot the 'raw' force
    ax1 = plt.subplot(n_rows,n_cols,1)
    plt.plot(time_plot,force_plot,**raw_force_kwargs)
    plt.plot(time_plot,force_interp_plot,**interp_force_kwargs)
    PlotUtilities.lazyLabel("","Force (pN)","")
    plt.ylim(ylim_force_pN)
    # plot arrows above the events
    arrow_kwargs = dict(plot_x=time_plot,plot_y=force_plot,
                        markersize=50)
    Plotting.plot_arrows_above_events(event_idx=info_final.event_idx,
                                      fudge_y=20,**arrow_kwargs)
    # add a scale bar
    scale_fraction_width = 0.3
    scale_fraction_offset = 0.3
    max_time = max(time_plot)
    width = scale_fraction_width * max_time
    label = "{:.1g}s".format(width)
    PlotUtilities.scale_bar_x(x=scale_fraction_offset*max_time,
                              y=-15,s=label,
                              width=width)     
    PlotUtilities.no_x_anything(ax1)               
    # # plot the 'zoomed' axis
    ax_zoom = plt.subplot(n_rows,n_cols,3)
    # determine the second event (zoom index)
    zoom_event_idx = 1
    event_zoom = info_final.event_idx[zoom_event_idx]
    zoom_slice = info_final.event_slices[zoom_event_idx]
    # get the zoom slice
    n_around_idx = 3e-3
    n_around = int(np.ceil(force_plot.size*n_around_idx))
    zoom_slice = slice(event_zoom-n_around,event_zoom+n_around)
    time_slice = time_plot[zoom_slice]
    force_slice = force_plot[zoom_slice]
    # get the fit to before the event (symmetric about the event location.)
    slice_fit = slice(0,int(time_slice.size/2),1)
    x_fit = time_slice[slice_fit]
    y_fit = force_slice[slice_fit]
    _,predicted,loading_rate,rupture_force,_ = \
        Analysis._loading_rate_helper(x_fit,y_fit)
    plt.plot(time_slice,force_slice,**raw_force_kwargs)
    color_loading_line = 'm'
    plt.plot(x_fit,predicted,linestyle='--',color=color_loading_line,linewidth=3)
    xlim = [time_slice[0],time_slice[-1]]
    plt.xlim(xlim)
    PlotUtilities.zoom_effect01(ax1, ax_zoom, *xlim)
    # plot annotations showing the loading rate
    bbox_props = dict(boxstyle="rarrow,pad=0.3",linestyle='--',
                      fc=color_loading_line,ec=color_loading_line,
                      alpha=0.3,lw=2)
    dy = (predicted[-1]-predicted[0])
    dx = (x_fit[-1]-x_fit[0])
    # XXX should programmtically figure out rotation...
    rotation =25
    fontsize=17
    t = ax_zoom.text(np.mean(x_fit), np.mean(force_slice), 
                     r"$\frac{dF}{dt}$", ha="center", va="center", 
                     rotation=rotation,size=fontsize,bbox=bbox_props)
    # plot annotations showing the rupture force
    fudge = (max(x_fit)-min(x_fit)) * 0.2
    min_time = min(x_fit) + fudge
    event_force = rupture_force
    event_time = time_plot[event_zoom]
    ax_zoom.annotate(r"F$_r$",
                     xy=(event_time,event_force), xycoords='data',
                     xytext=(min_time, event_force), textcoords='data',
                     verticalalignment='center',
                     horizontalalignment='center',fontsize=fontsize,
                     arrowprops=dict(arrowstyle="->",
                                     connectionstyle="arc3"))
    zoom_event_only = [event_zoom]
    Plotting.plot_arrows_above_events(event_idx=zoom_event_only,fudge_y=6,
                                      **arrow_kwargs)
    # add a scalebar...
    dx_zoom_full =abs(time_slice[-1]-time_slice[0])
    width = scale_fraction_width * dx_zoom_full
    label = "{:.1g}ms".format(1000*width)
    x_text = time_slice[0] + dx_zoom_full*scale_fraction_offset
    PlotUtilities.scale_bar_x(x=x_text,
                              y=min(force_slice)*1.1,s=label,
                              width=width)    
    PlotUtilities.lazyLabel("Time","Force (pN)","")    
    PlotUtilities.no_x_anything(ax_zoom)                                   
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
