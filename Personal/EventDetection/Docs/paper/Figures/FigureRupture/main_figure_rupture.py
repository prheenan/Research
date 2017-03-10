# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection.Util import Analysis 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

def fmt(ax):
    ax.set_xlim([-0.1,2.5])
    ax.set_ylim([-50,40])


def before_and_after(x,y,before_slice,after_slice,style):
    color_before = 'b'
    color_after = 'r'
    tuples = [ [x,y,before_slice,color_before,style],
               [x,y,after_slice,color_after,style]]
    for x_tmp,y_tmp,slice_v,color_tmp,style_tmp in tuples:
        x_sliced = x_tmp[slice_v]
        plt.plot(x_sliced,y_tmp[slice_v],color=color_tmp,**style_tmp)


def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    out_fig = "out.png"
    example = read_and_cache_file(data_base + "rupture.csv",has_events=True,
                                  force=False,cache_directory=data_base)
    n_filter = 1000
    kw = dict(cache_directory=data_base,force=False)
    events = example.Events
    fec_split = Analysis.zero_and_split_force_extension_curve(example)
    event_idx_retract = fec_split.get_retract_event_centers()
    points_around = 2000
    retract = fec_split.retract
    retract.Force -= np.median(retract.Force)
    retract_filtered = FEC_Util.GetFilteredForce(retract,n_filter)
    # get everything in terms of ploting variables
    x_plot = lambda x: x - min(x)
    y_plot = lambda y: y * 1e12
    x = x_plot(retract.Time)
    force = y_plot(retract.Force)
    x_filtered = x_plot(retract_filtered.Time)
    force_filtered = y_plot(retract_filtered.Force)
    event_idx = event_idx_retract[0]
    zoom_start_idx = event_idx-points_around
    zoom_end_idx = event_idx+points_around
    slice_event_subplot = slice(zoom_start_idx,zoom_end_idx,1)
    # determine where the rupture actually occured
    slice_before_zoom = slice(zoom_start_idx,event_idx)
    slice_after_zoom = slice(event_idx,zoom_end_idx)
    x_event = x[slice_before_zoom]
    f_event = force[slice_before_zoom]
    coeffs,predicted,loading_rate,rupture_force,index = \
        Analysis._loading_rate_helper(x_event,f_event)
    index_absolute = index + slice_event_subplot.start
    rupture_time = x[index_absolute]
    # update all the slices
    slice_before_zoom = slice(zoom_start_idx,index_absolute)
    slice_after_zoom = slice(index_absolute,zoom_end_idx)
    slice_before_event = slice(0,index_absolute,1)
    slice_after_event = slice(index_absolute,None,1)
    x_zoom = x[slice_event_subplot]
    f_zoom = force[slice_event_subplot]
    xlim = [min(x_zoom),max(x_zoom)]
    ylim = [min(f_zoom),max(f_zoom)]
    fig = PlotUtilities.figure((16,8))
    plt.subplot(1,2,1)
    ax = plt.gca()
    fmt(ax)
    style_data = dict(alpha=0.3,linewidth=1)
    style_filtered = dict(alpha=1,linewidth=2)
    # plot the force etc
    before_and_after(x,force,slice_before_event,slice_after_event,style_data)
    before_and_after(x,force_filtered,slice_before_event,slice_after_event,
                     style_filtered)
    # add a rectangle and its border
    min_f = min(force)
    max_f = max(force)
    plt.ylim([min_f,max_f])
    norm = lambda x: (x-min_f)/(max_f-min_f)
    ymin_box = norm(min(ylim))
    ymax_box = norm(max(ylim))
    y_args = dict(ymin=ymin_box,ymax=ymax_box)
    plt.axvspan(*xlim, fill=False,linestyle='dashdot',edgecolor='k',
                linewidth=3,**y_args)
    plt.axvspan(*xlim, fill=True,color='k',alpha=0.15,
                **y_args)
    PlotUtilities.lazyLabel("Time [s]","Force (pN)","")
    # plot the rupture
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    # zoom-factor: 2.5, location: upper-left
    plt.subplot(1,2,2)
    before_and_after(x,force,slice_before_zoom,slice_after_zoom,style_data)
    plt.plot(x_event,predicted,color='k',linestyle='--',linewidth=3)
    plt.plot(x[index_absolute],predicted[index],'go',markersize=10,linewidth=0,
             alpha=0.3)
    plt.xlim(xlim)
    plt.ylim(ylim)
    PlotUtilities.lazyLabel("Time [s]","","")
    PlotUtilities.savefig(fig,out_fig)
    

if __name__ == "__main__":
    run()
