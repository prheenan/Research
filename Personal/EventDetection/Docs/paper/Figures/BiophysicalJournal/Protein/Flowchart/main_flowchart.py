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
    force_plot = split_fec_no_domain_specific.retract.Force * 1e12
    force_interp_plot = \
        split_fec_no_domain_specific.retract_spline_interpolator()(time) * 1e12
    # plot everything
    raw_force_kwargs = dict(color='k',alpha=0.3)
    interp_force_kwargs = dict(color='b',linewidth=1.5)
    probabiity_kwargs = dict(color='r',linestyle=":")
    n_cols = 1
    n_rows = 6
    ylim_force_pN = [-30,max(force_interp_plot)*1.1]
    to_prob_plot = lambda x: np.log10(x)
    print(info_final.cdf,to_prob_plot(info_final.cdf))
    ylim_prob = [to_prob_plot(min((info_final.cdf))/5),1.1]
    title_kwargs = dict(loc='left')
    kw = dict(title_kwargs=title_kwargs)
    arrow = "$\downarrow$"
    probability_label = "$\log_{10}$(Probability)"
    fig = PlotUtilities.figure(figsize=(8,16))
    # plot the 'raw' force
    ax_raw = plt.subplot(n_rows,n_cols,1)
    plt.plot(time_plot,force_plot,**raw_force_kwargs)
    PlotUtilities.lazyLabel("","Force (pN)","")
    plt.ylim(ylim_force_pN)
    PlotUtilities.no_x_label(ax_raw)
    # plot the spline fit
    ax_spline = plt.subplot(n_rows,n_cols,2)
    plt.plot(time_plot,force_plot,**raw_force_kwargs)    
    plt.plot(time_plot,force_interp_plot,**interp_force_kwargs)
    title_spline = arrow + " Spline fit"
    PlotUtilities.lazyLabel("","Force (pN)",title_spline,**kw)
    plt.ylim(ylim_force_pN)                            
    PlotUtilities.no_x_label(ax_spline)    
    # plot the 'raw' probability
    ax_raw_prob = plt.subplot(n_rows,n_cols,3)
    plt.plot(time_plot,to_prob_plot(info_no_domain_specific.cdf),
             **probabiity_kwargs)    
    title_prob = arrow + "Determine probabilty of non-event"
    PlotUtilities.lazyLabel("",probability_label,title_prob,**kw)
    PlotUtilities.no_x_label(ax_raw_prob)        
    plt.ylim(ylim_prob)
    # plot the adhesion-fixed probability
    ax_adhesion = plt.subplot(n_rows,n_cols,4)
    plt.plot(time_plot,to_prob_plot(info_remove_adhesions.cdf),
             **probabiity_kwargs)    
    title_adhesion = arrow + "Supress force increases and surface adhesion"
    PlotUtilities.lazyLabel("","",title_adhesion,**kw)
    PlotUtilities.no_x_label(ax_adhesion)      
    plt.ylim(ylim_prob)    
    # plot the final probability
    ax_final_prob = plt.subplot(n_rows,n_cols,5)
    plt.plot(time_plot,to_prob_plot(info_final.cdf),
             **probabiity_kwargs)    
    title_consistent = (arrow + "Supress where force change consistent with 0")
    PlotUtilities.lazyLabel("","",title_consistent,**kw)
    PlotUtilities.no_x_label(ax_final_prob)      
    plt.ylim(ylim_prob)    
    # plot the final event locations
    ax = plt.subplot(n_rows,n_cols,6)
    plt.plot(time_plot,force_plot,**raw_force_kwargs)    
    plt.plot(time_plot,force_interp_plot,**interp_force_kwargs)
    title_final = (arrow + " Extract events where probability < threshold")
    PlotUtilities.lazyLabel("Time (s)","Force (pN)",title_final,**kw)
    for e in info_final.event_slices_raw:
        time_slice = time_plot[e]
        bounds = [time_slice[0],time_slice[-1]]
        plt.axvspan(*bounds,alpha=0.3,color='r')
    plt.ylim(ylim_force_pN)                           
    event_idx_fudge_and_kw = \
        [ [0 ,0.005,dict(loc=2,inset_kwargs=dict(loc1=3,loc2=4))],
          [-1,0.007,dict(loc=1,inset_kwargs=dict(loc1=2,loc2=3))] ]
    for i,fudge,kw in event_idx_fudge_and_kw:
        # determine the slice we want to use                 
        event_location = info_final.event_idx[i]
        event_slice = info_final.event_slices_raw[i]
        extra_for_inset_plot = -int(np.ceil(time_plot.size * fudge))
        event_bounding_slice = info_final.event_slices[i]
        # add in some padding 
        event_bounding_slice = \
            slice(event_bounding_slice.start-extra_for_inset_plot,
                  event_bounding_slice.stop+extra_for_inset_plot,1)
        time_first_event_plot = time_plot[event_bounding_slice]       
        inset_by_slice(ax,time_plot,force_plot,force_interp_plot,
                       event_bounding_slice,
                       event_location,raw_force_kwargs,
                       interp_force_kwargs,**kw)
    # set all the labels how we want
    axs = fig.axes
    for ax_tmp in axs:
        yax = ax_tmp.get_yaxis()
        yax.set_label_coords(-0.1,0.5)
    PlotUtilities.savefig(fig,"./flowchart.png")
    
def inset_by_slice(ax,x,y,y_interp,slice_v,event_idx,y_kwargs,interp_kwargs,loc,
                   inset_kwargs=dict()):
    # add an inset box for the first (assumed hardest-to-find) event
    in_ax = inset_axes(ax,
                       width="30%", # width = X% of parent_bbox
                       height="50%", # height  = X% of parent_bbox
                       loc=loc)    
   # plot the raw and interpolated force
    mark_inset(ax, in_ax, fc="none", ec="0.5",**inset_kwargs)    
    in_ax.axvline(x[event_idx])    
    in_ax.plot(x[slice_v],y[slice_v],**y_kwargs)
    in_ax.plot(x[slice_v],y_interp[slice_v],**interp_kwargs)       
    PlotUtilities.no_x_anything(ax=in_ax)
    PlotUtilities.no_y_anything(ax=in_ax)                       
    

if __name__ == "__main__":
    run()
