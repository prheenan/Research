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
from Research.Personal.EventDetection.Util import Analysis,Plotting
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
    BboxConnectorPatch
from matplotlib.transforms import Bbox, TransformedBbox, \
    blended_transform_factory



def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)
    bbox_patch1 = BboxPatch(bbox1, color='k',**prop_patches)
    bbox_patch2 = BboxPatch(bbox2, color='w',**prop_patches)
    p = BboxConnectorPatch(bbox1, bbox2,
                           # loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes
    (xmin,xmax) : the limits of the colored area in both plot axes.

    connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
    be marked.  The keywords parameters will be used ti create
    patches.

    """

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches = kwargs.copy()
    alpha = 0.2
    prop_patches["ec"] = "none"
    prop_patches["alpha"] = alpha
    prop_lines = dict(color='k',alpha=alpha,**kwargs)
    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=prop_lines, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    out_fig = "ruptures.pdf"
    example = read_and_cache_file(data_base + "rupture.csv",has_events=True,
                                  force=False,cache_directory=data_base)
    n_filter = 1000
    kw = dict(cache_directory=data_base,force=False)
    events = example.Events
    fec_split = Analysis.zero_and_split_force_extension_curve(example)
    event_idx_retract = fec_split.get_retract_event_centers()
    event_idx = event_idx_retract[0]
    zoom_factor = 20
    points_around = int(np.ceil(event_idx/zoom_factor))
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
    zoom_start_idx = event_idx-points_around
    zoom_end_idx = event_idx+points_around
    slice_event_subplot = slice(zoom_start_idx,zoom_end_idx,1)
    # determine where the rupture actually occured
    slice_before_zoom = slice(zoom_start_idx,event_idx)
    slice_after_zoom = slice(event_idx,zoom_end_idx)
    x_event = x[slice_before_zoom]
    f_event = force[slice_before_zoom]
    x_event_both = x[slice_event_subplot]
    f_event_both = force[slice_event_subplot]
    coeffs,predicted,loading_rate,rupture_force,index = \
        Analysis._loading_rate_helper(x_event,f_event)
    predicted_both = np.polyval(coeffs,x_event_both)
    index_after_event = index+1
    index_absolute = index_after_event + slice_event_subplot.start
    rupture_time = x[index_absolute]
    # update all the slices
    slice_before_zoom = slice(zoom_start_idx,index_absolute)
    slice_after_zoom = slice(index_absolute,zoom_end_idx)
    slice_before_event = slice(0,index_absolute,1)
    slice_after_event = slice(index_absolute,None,1)
    x_zoom = x[slice_event_subplot]
    f_zoom = force[slice_event_subplot]
    ylim = [-30,max(force)*2]
    xlim = [0,1.0]
    xlim_zoom = [min(x_zoom),max(x_zoom)]
    ylim_zoom = [min(f_zoom),max(f_zoom)]
    # get the second zoom
    second_zoom = int(points_around/zoom_factor)
    zoom_second_before = slice(index_absolute-second_zoom,index_absolute)
    zoom_second_after  = slice(index_absolute,index_absolute+second_zoom)
    slice_second_zoom = slice(zoom_second_before.start,
                              zoom_second_after.stop)
    x_second_zoom = x[slice_second_zoom]
    force_second_zoom = force[slice_second_zoom]
    xlim_second_zoom = [min(x_second_zoom),max(x_second_zoom)]
    ylim_second_zoom = [min(force_second_zoom),max(force_second_zoom)]
    # plot everything
    n_rows = 3
    n_cols = 1
    fudge_y = 5
    ylabel = "Force (pN)"
    fig = PlotUtilities.figure((8,12))
    ax1 = plt.subplot(n_rows,n_cols,1)
    style_data = dict(alpha=0.5,linewidth=1)
    style_filtered = dict(alpha=1.0,linewidth=3)
    # plot the force etc
    Plotting.before_and_after(x,force,slice_before_event,slice_after_event,
                              style_data,label="Raw (25kHz)")
    Plotting.before_and_after(x,force_filtered,slice_before_event,
                              slice_after_event,style_filtered,
                              label="Filtered (25Hz)")
    PlotUtilities.lazyLabel("Time (s)",ylabel,"",loc="upper left",
                            frameon=False)
    plt.ylim(ylim)
    plt.xlim(xlim)
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top') 
    Plotting.plot_arrows_above_events([index_absolute],x,force_filtered,
                                      fudge_y=fudge_y*4.5)
    # plot the rupture
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    # zoom-factor: 2.5, location: upper-left
    ax2 = plt.subplot(n_rows,n_cols,2)
    style_data_post_zoom = dict(style_data)
    post_alpha = style_data['alpha']
    style_data_post_zoom['alpha']=post_alpha
    style_filtered_post_zoom = dict(style_filtered)
    style_filtered_post_zoom['alpha']=1
    Plotting.before_and_after(x,force,slice_before_zoom,slice_after_zoom,
                              style_data_post_zoom)
    plot_rupture = lambda l: Plotting.\
        plot_arrows_above_events([index_after_event],x_event_both,
                                 predicted_both,fudge_y=fudge_y)
    plot_line = lambda l :  plt.plot(x_event_both[:index_after_event+1],
                                     predicted_both[:index_after_event+1],
                                     color='m',
                                     linestyle='-',linewidth=3,label=l)
    plot_rupture('Rupture')
    plot_line("Fit")
    PlotUtilities.lazyLabel("",ylabel,"",frameon=False,
                            loc='upper right',legend_kwargs=dict(numpoints=1))
    plt.xlim(xlim_zoom)
    plt.ylim(ylim_zoom)
    # make a scale bar for this plot
    scale_width = 0.01
    string = "{:d} ms".format(int(scale_width*1000))
    get_bar_location = lambda _xlim: np.mean([np.mean(_xlim),min(_xlim)])
    PlotUtilities.scale_bar_x(get_bar_location(xlim_zoom),
                              -10,s=string,width=scale_width)
    PlotUtilities.no_x_label()
    ax3 = plt.subplot(n_rows,n_cols,3)
    Plotting.before_and_after(x,force,zoom_second_before,zoom_second_after,
                              style_data_post_zoom)
    plt.xlim(xlim_second_zoom)
    plt.ylim(ylim_second_zoom)
    plot_rupture("")
    plot_line("")
    PlotUtilities.lazyLabel("",ylabel,"")
    # make a scale bar for this plot
    scale_width = scale_width/zoom_factor
    string = "{:.1f} ms".format(scale_width*1000)
    PlotUtilities.scale_bar_x(get_bar_location(xlim_second_zoom),
                              0,s=string,width=scale_width)
    PlotUtilities.no_x_label()
    # draw lines connecting the plots
    zoom_effect01(ax1, ax2, *xlim_zoom,linewidth=3)
    zoom_effect01(ax2, ax3, *xlim_second_zoom,linewidth=3)
    subplots_adjust = dict(hspace=0.1)
    # save out without the labels
    PlotUtilities.savefig(fig,out_fig.replace(".pdf","_pres.pdf"),
                          subplots_adjust=subplots_adjust,close=False)
    # save out with the labels
    PlotUtilities.label_tom(fig,loc=(-0.1,0.97))
    PlotUtilities.savefig(fig,out_fig,subplots_adjust=subplots_adjust)
    

if __name__ == "__main__":
    run()
