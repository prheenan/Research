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
import matplotlib.gridspec as gridspec

def make_scale_bar(y_frac=0.2,x_frac=0.2,width=0.2,unit="s",mult=1):
    xlim = plt.xlim()
    x_full = abs(np.diff(xlim))[0]
    width = width * x_full
    ylim = plt.ylim()
    y = np.max(ylim) - abs(np.diff(ylim))*y_frac
    x = np.min(xlim) + abs(x_full) * x_frac
    s = "{:.1g}{:s}".format(width*mult,unit)
    PlotUtilities.scale_bar_x(x=x,y=y,s=s,width=width,fontsize=18)

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
    # get the 'raw' no-event probabilities, and the increasingly domain-specific
    # ones
    threshold = 1e-3
    pred_kw = dict(threshold=threshold)
    split_fec_no_domain_specific,info_no_domain_specific= \
        Detector._predict_full(data,f_refs=[],**pred_kw)
    _,info_remove_adhesions = Detector._predict_full(data,\
      f_refs=[Detector.adhesion_mask_function_for_split_fec],**pred_kw)
    _,info_final = Detector._predict_full(data,**pred_kw)
    # get the plotting versions of the time, etc
    time= split_fec_no_domain_specific.retract.Time
    time_plot = time - time[0]
    final_event_time = time_plot[info_final.event_idx[-1]]
    xlim_time = [0,final_event_time*1.1]
    force_plot = split_fec_no_domain_specific.retract.Force * 1e12
    force_interp_plot = \
        split_fec_no_domain_specific.retract_spline_interpolator()(time) * 1e12
    # plot everything
    raw_force_kwargs = dict(color='k',alpha=0.3)
    interp_force_kwargs = dict(color='b',linewidth=1.5)
    probabiity_kwargs = dict(color='r',linestyle="-")
    # for the probabilities, how far from the maximum y the scale bar should be
    # (units [0,1]
    y_frac_prob = 0.7
    n_cols = 1
    n_rows = 6
    ylim_force_pN = [-30,max(force_interp_plot)+1.1+35]
    to_prob_plot = lambda x: np.log10(x)
    ylim_prob = [to_prob_plot(min((info_final.cdf))/5),1.1]
    title_kwargs = dict(loc='left',fontsize=18,color='b')
    kw = dict(title_kwargs=title_kwargs)
    arrow = "$\downarrow$"
    probability_label = "$\log_{10}$($P$)"
    probability_label_post = ""
    n_cols = 3
    n_rows = 6
    gs = gridspec.GridSpec(nrows=n_rows,ncols=n_cols,
                           width_ratios=[1 for _ in range(n_cols)],
                           height_ratios=[2,1,1,1,3,3])
    fig = PlotUtilities.figure(figsize=(5,10))
    # plot the 'raw' force and spline
    ax_raw = plt.subplot(gs[0,:])
    plt.plot(time_plot,force_plot,label="Raw",**raw_force_kwargs)    
    plt.plot(time_plot,force_interp_plot,label="Spline",
             **interp_force_kwargs)
    PlotUtilities.x_label_on_top(ax_raw)
    PlotUtilities.no_x_label(ax_raw)
    plt.ylim(ylim_force_pN)                
    PlotUtilities.lazyLabel("Time (s)","Force (pN)","",loc="upper center",
                            legend_kwargs=dict(handlelength=1,fontsize=16))
    plt.xlim(xlim_time)
    make_scale_bar()
    # # plot the 'raw' probability
    ax_raw_prob = plt.subplot(gs[1,:])
    plt.plot(time_plot,to_prob_plot(info_no_domain_specific.cdf),
             **probabiity_kwargs)    
    title_prob = arrow + "Get probability of no event"
    PlotUtilities.lazyLabel("",probability_label,title_prob,**kw)
    PlotUtilities.no_x_label(ax_raw_prob)        
    plt.ylim(ylim_prob)
    plt.xlim(xlim_time)
    make_scale_bar(y_frac=y_frac_prob)
    # # plot the adhesion-fixed probability
    ax_adhesion = plt.subplot(gs[2,:])
    plt.plot(time_plot,to_prob_plot(info_remove_adhesions.cdf),
             **probabiity_kwargs)    
    title_adhesion = arrow + r"Supress adhesion and $\frac{dF}{dt}>0$"
    PlotUtilities.lazyLabel("",probability_label_post,title_adhesion,**kw)
    PlotUtilities.no_x_label(ax_adhesion)      
    plt.ylim(ylim_prob)
    plt.xlim(xlim_time)    
    make_scale_bar(y_frac=y_frac_prob)
    # # plot the final probability
    ax_final_prob = plt.subplot(gs[3,:])
    plt.plot(time_plot,to_prob_plot(info_final.cdf),
             **probabiity_kwargs)    
    title_consistent = (arrow + "Supress force changes near 0")
    PlotUtilities.lazyLabel("",probability_label_post,title_consistent,**kw)
    PlotUtilities.no_x_label(ax_final_prob)      
    plt.ylim(ylim_prob)    
    plt.xlim(xlim_time)
    make_scale_bar(y_frac=y_frac_prob)
    # # plot the final event locations
    ax_final = plt.subplot(gs[4,:])
    plt.plot(time_plot,force_plot,**raw_force_kwargs)    
    plt.plot(time_plot,force_interp_plot,**interp_force_kwargs)
    PlotUtilities.no_x_label(ax_final)      
    title_final = (arrow + " Extract significant events")
    event_starts = [e.start for e in info_final.event_slices_raw]
    Plotting.plot_arrows_above_events(event_starts,plot_x=time_plot,
                                      plot_y=force_plot,fudge_y=25,
                                      label=None)
    PlotUtilities.lazyLabel("","Force (pN)",title_final,
                            loc = "upper center",**kw)
    plt.ylim(ylim_force_pN)                           
    plt.xlim(xlim_time)
    make_scale_bar()
    ylim_first_event = [-5,30]
    first_event_window_large = 0.02
    fraction_increase= 5
    # get the event index, window pct to use, where to show a 'zoom', and the
    # y limits (if none, just all of it)
    event_idx_fudge_and_kw = \
        [ [0 ,first_event_window_large   ,True,ylim_first_event],
          [0 ,first_event_window_large/fraction_increase,False,
           ylim_first_event],
          [-1,4e-3,True,[-50,None]]]
    for i,(event_id,fudge,zoom_bool,ylim) in enumerate(event_idx_fudge_and_kw):
        # determine the slice we want to use                 
        event_location = info_final.event_idx[event_id]
        event_slice = info_final.event_slices_raw[event_id]
        extra_for_inset_plot = int(np.ceil(time_plot.size * fudge))
        event_bounding_slice = info_final.event_slices[event_id]
        # add in some padding 
        event_bounding_slice = \
            slice(event_location-extra_for_inset_plot,
                  event_location+extra_for_inset_plot,1)
        time_first_event_plot = time_plot[event_bounding_slice]
        in_ax = plt.subplot(gs[-1,i])
        time_slice = time_first_event_plot
        in_ax.plot(time_slice,force_plot[event_bounding_slice],
                   **raw_force_kwargs)
        in_ax.plot(time_slice,force_interp_plot[event_bounding_slice],
                   **interp_force_kwargs)       
        PlotUtilities.no_x_anything(ax=in_ax)
        if (i == 0):
            ylabel = "Force (pN)"
        else:
            ylabel = ""
        PlotUtilities.lazyLabel("Time (s)",ylabel,"")
        # determine if we need to add in 'guidelines' for zooming
        if (zoom_bool):
            PlotUtilities.zoom_effect01(ax_final,in_ax,*in_ax.get_xlim())
        else:
            # this is a 'second' zoom in...'
            PlotUtilities.no_y_label(in_ax)
            ylab = (r"$\downarrow${:d}x").format(fraction_increase)
            PlotUtilities.ylabel(ylab)
        # plot an arrow over the (single) event
        Plotting.plot_arrows_above_events([event_location],plot_x=time_plot,
                                          plot_y=force_plot,fudge_y=7,
                                          label=None,markersize=150)
        plt.ylim(ylim)
        make_scale_bar(y_frac=0.9,x_frac=0.5,width=0.2,unit="ms",mult=1000)
    PlotUtilities.savefig(fig,"./flowchart.png",
                          subplots_adjust=dict(hspace=0.35))
    

if __name__ == "__main__":
    run()
