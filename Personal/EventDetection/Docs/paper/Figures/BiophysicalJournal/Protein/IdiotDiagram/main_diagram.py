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

import scipy

sys.path.append("../../../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting,Learning,\
    Analysis
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from mpl_toolkits.axes_grid.inset_locator import inset_axes,mark_inset

def escape_rate(loading_rate,rupture_force,delta_G_dagger,x_dagger,k0,
                nu,beta):
    a = (1-nu*rupture_force*x_dagger/delta_G_dagger)
    b = 1-a**(1/nu)
    return k0 * (a ** (1/nu -1)) * np.exp(beta*delta_G_dagger * b)

def dudko_model(loading_rate,rupture_force,delta_G_dagger,x_dagger,k0,nu,beta):
    k_F = escape_rate(loading_rate,rupture_force,delta_G_dagger,x_dagger,
                      k0,nu,beta)
    c = beta*x_dagger*loading_rate                      
    d = (1-nu * rupture_force * x_dagger/delta_G_dagger)**(1-1/nu)
    return (1/loading_rate) * k_F * np.exp(k0/c) * np.exp( (-k_F/c) * d)
            
def mean_rupture_force(loading_rate,delta_G_dagger,x_dagger,k0,nu,beta):
    c0 = delta_G_dagger*beta
    gamma = 0.577
    f = (1/c0) * np.log(k0*np.exp(c0+gamma)/(beta*x_dagger*loading_rate))
    return delta_G_dagger/(nu*x_dagger) * (1- f**nu)
    
def stdev_rupture_force(loading_rate,delta_G_dagger,x_dagger,k0,nu,beta):
    c0 = delta_G_dagger*beta
    gamma_t = 1.064
    f = np.log(k0*np.exp(c0+gamma_t)/(beta*x_dagger*loading_rate))
    variance = (np.pi**2/(6*(beta*x_dagger)**2)) * (1/c0 * f)**(2*nu-2)
    return np.sqrt(variance)
            
def  generate_normalized_model(loading_rate,rupture_forces,**kwargs):
    model = dudko_model(loading_rate,rupture_forces,**kwargs)
    model = model/sum(model)                            
    return model
    
def generate_rupture_histogram():
    np.random.seed(42)
    kbT =  4.1e-21
    delta_G_dagger = 10 *kbT
    x_dagger= 0.3e-9
    k0 = 0.1
    nu = 2/3
    beta = 1/kbT
    n_samples = 1000
    # loading rates from 10 pN/s to 100 pN/s, 
    loading_rates = np.logspace(-11,-10,num=5)
    # rupture forces from 1pN to 1000pN
    F_c = delta_G_dagger/(nu*x_dagger)
    rupture_forces = np.linspace(1e-12,F_c,num=500)
    common_kwargs = dict(delta_G_dagger=delta_G_dagger,
                         x_dagger=x_dagger,k0=k0,beta=beta,
                         nu=nu)
    mean_rupture_forces = mean_rupture_force(loading_rates,**common_kwargs)
    stdev_rupture_forces = stdev_rupture_force(loading_rates,**common_kwargs)
    models = []    
    for loading_rate in loading_rates:
        model = generate_normalized_model(loading_rate,rupture_forces,
                                          **common_kwargs)
        models.append(model)                                        
    # generate histograms from all the models 
    rupture_forces_histograms = []
    for loading_rate,model in zip(loading_rates,models):
        # ensure they sum to 1 (XXX shouldnt be needed?)
        choices = np.random.choice(a=rupture_forces,size=n_samples,p=model)
        rupture_forces_histograms.append(choices)
    return loading_rates,rupture_forces_histograms,models,\
        rupture_forces,mean_rupture_forces,stdev_rupture_forces

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    loading_rate_histogram,rupture_forces_histograms,models,rupture_forces,\
        mean_rupture_forces,stdev_rupture_forces = generate_rupture_histogram() 
    loading_plot = loading_rate_histogram*1e12
    for rupture,loading in zip(rupture_forces_histograms,
                               loading_plot):
        rupture_plot = rupture * 1e12                               
        most_likely_rupture = np.mean(rupture_plot)
        std_rupture = np.std(rupture_plot)
        plt.errorbar(loading,y=most_likely_rupture,yerr=std_rupture,fmt='ro')
        PlotUtilities.lazyLabel("Loading rate(pN/s)",
                                "Mean Rupture Force(pN)","")                                  
    mean_plot = mean_rupture_forces*1e12                                
    stdev_plot = stdev_rupture_forces*1e12
    plt.plot(loading_plot,mean_plot)
    plt.plot(loading_plot,mean_plot-stdev_plot,linestyle='--')
    plt.plot(loading_plot,mean_plot+stdev_plot,linestyle='--')    
    plt.show()
    for rupture,model in zip(rupture_forces_histograms,models):
        n,_,_ = plt.hist(rupture*1e12)
        plt.plot(rupture_forces*1e12,model*max(n)/max(model))
        PlotUtilities.lazyLabel("Force (pN)","Count","")
        plt.show()
    exit(1)
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
    PlotUtilities.savefig(fig,"./diagram.png")

if __name__ == "__main__":
    run()
