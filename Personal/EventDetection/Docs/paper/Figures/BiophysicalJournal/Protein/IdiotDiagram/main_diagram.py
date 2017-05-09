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
from matplotlib.patches import Ellipse
from Research.Personal.EventDetection.Util import Offline,Plotting,Learning,\
    Analysis
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from mpl_toolkits.axes_grid.inset_locator import inset_axes,mark_inset

from FitUtil.EnergyLandscapes.Rupture_Dudko2007.Python.Code import Dudko2007
import matplotlib.gridspec as gridspec

# plotting constants    
raw_force_kwargs = dict(color='k',alpha=0.3)
interp_force_kwargs = dict(color='b',linewidth=3)
probabiity_kwargs = dict(color='r')
# how big the scale bars are
scale_fraction_width = 0.3
scale_fraction_offset = 0.3
df_dt_string = r"$\frac{\mathrm{dF}}{\mathrm{dt}}$"
rupture_string = r"F$_{\mathrm{r}}$"
fontsize=20

def generate_rupture_histograms():
    num_loading_rates = 10
    np.random.seed(42)
    kbT =  4.1e-21
    delta_G_ddagger = 10 *kbT
    x_ddagger= 0.3e-9
    k0 = 0.1
    nu = 2/3
    beta = 1/kbT
    n_samples = 10000
    # loading rates from 10 pN/s to 100 pN/s, 
    loading_rates = np.logspace(-11,-10,num=num_loading_rates)
    # rupture forces from 1pN to 1000pN
    F_c = delta_G_ddagger/(nu*x_ddagger)
    rupture_forces = np.linspace(1e-12,F_c,num=n_samples)
    basic_kwargs = dict(delta_G_ddagger=delta_G_ddagger,
                        x_ddagger=x_ddagger,nu=nu)
    common_kwargs = dict(k0=k0,beta=beta,**basic_kwargs)
    mean_rupture_forces = Dudko2007.mean_rupture_force(loading_rates,
                                                       **common_kwargs)
    stdev_rupture_forces = Dudko2007.stdev_rupture_force(loading_rates,
                                                         **common_kwargs)
    models = []    
    x = np.linspace(start=-x_ddagger*0.8,stop=x_ddagger*0.8)
    landscape = Dudko2007.free_energy_landscape(x=x,**basic_kwargs)
    for loading_rate in loading_rates:
        model = Dudko2007.normalized_model(loading_rate,rupture_forces,
                                           **common_kwargs)
        models.append(model)                                        
    # generate histograms from all the models 
    rupture_forces_histograms = []
    for loading_rate,model in zip(loading_rates,models):
        # ensure they sum to 1 (XXX shouldnt be needed?)
        choices = np.random.choice(a=rupture_forces,size=n_samples,p=model)
        rupture_forces_histograms.append(choices)
    return loading_rates,rupture_forces_histograms,models,\
        rupture_forces,mean_rupture_forces,stdev_rupture_forces,x,landscape

def plot_fec_scaled(time_plot,force_plot,force_interp_plot,info_final,
                    arrow_kwargs):
    plt.plot(time_plot,force_plot,**raw_force_kwargs)
    plt.plot(time_plot,force_interp_plot,**interp_force_kwargs)
    PlotUtilities.lazyLabel("","Force (pN)","")
    # plot arrows above the events
    Plotting.plot_arrows_above_events(event_idx=info_final.event_idx,
                                      fudge_y=20,**arrow_kwargs)
    # add a scale bar
    PlotUtilities.no_x_anything(plt.gca())               
    max_time = max(time_plot)
    width = scale_fraction_width * max_time
    label = "{:.1g}s".format(width)
    PlotUtilities.scale_bar_x(x=scale_fraction_offset*max_time,
                              y=-15,s=label,
                              width=width)  

def common_arrow_kwargs(arrowprops=dict(arrowstyle="<->",shrinkA=20,
                                        shrinkB=20,
                                        connectionstyle="arc3")):
    return  dict(xycoords='data',
                 textcoords='data',
                 verticalalignment='center',
                 horizontalalignment='center',fontsize=fontsize,
                 arrowprops=arrowprops)

def plot_zoomed(time_plot,force_plot,info_final,ax1,arrow_kwargs):
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
    plt.plot(x_fit,predicted,linestyle='--',color=color_loading_line,
             linewidth=3)
    xlim = [time_slice[0],time_slice[-1]]
    plt.xlim(xlim)
    ax_zoom = plt.gca()
    PlotUtilities.zoom_effect01(ax1, ax_zoom, *xlim)
    # plot annotations showing the loading rate
    bbox_props = dict(boxstyle="rarrow,pad=0.3",linestyle='--',
                      fc=color_loading_line,ec=color_loading_line,
                      alpha=0.3,lw=2)
    dy = (predicted[-1]-predicted[0])
    dx = (x_fit[-1]-x_fit[0])
    # XXX should programmtically figure out rotation...
    rotation =25
    t = ax_zoom.text(np.mean(x_fit), np.mean(force_slice), 
                     df_dt_string, ha="center", va="center", 
                     rotation=rotation,size=fontsize,bbox=bbox_props)
    # plot annotations showing the rupture force
    fudge = (max(x_fit)-min(x_fit)) * 0.2
    min_time = min(x_fit) + fudge
    event_force = rupture_force
    event_time = time_plot[event_zoom]
    arrow_kw_thick = common_arrow_kwargs()
    arrow_kw_thick['arrowprops']['linewidth'] = 3
    ax_zoom.annotate(rupture_string,
                     xy=(event_time,event_force), 
                     xytext=(min_time, event_force),**arrow_kw_thick)
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


def plot_mean_rupture(rupture_forces_histograms,loading_rate_histogram,
                      mean_rupture_forces,stdev_rupture_forces):
    loading_plot = loading_rate_histogram*1e12
    means = np.mean(rupture_forces_histograms,axis=1) *1e12
    stdevs= np.std(rupture_forces_histograms,axis=1) * 1e12
    plt.errorbar(loading_plot,y=means,yerr=stdevs,fmt='ro')
    mean_plot = mean_rupture_forces*1e12                                
    stdev_plot = stdev_rupture_forces*1e12
    plt.plot(loading_plot,mean_plot)
    plt.plot(loading_plot,mean_plot-stdev_plot,linestyle=':')
    plt.plot(loading_plot,mean_plot+stdev_plot,linestyle=':')    

def plot_histogram_and_model(rupture_forces,rupture,model,kwargs_errorbar):
    rupture_plot = rupture*1e12
    n,_,_ = plt.hist(rupture_plot,alpha=0.5)
    plt.plot(rupture_forces*1e12,model*max(n)/max(model),linewidth=3)
    mean = np.mean(rupture_plot)
    stdev = np.std(rupture_plot)
    plt.errorbar(x=mean,xerr=stdev,y=max(n)*1.05,**kwargs_errorbar)

def plot_landscape(x,landscape):
    ax = plt.gca()
    landscape -= min(landscape)
    landscape /= 4.1e-21
    x *= 1e9
    plt.plot(x,landscape)
    PlotUtilities.no_x_label(ax)
    # determine where to put all the annotations
    max_idx = np.argmax(landscape)
    min_idx =np.argmin(landscape)
    x_max = x[max_idx]
    x_min = x[min_idx]
    y_low_plot = -2
    y_max = landscape[max_idx]
    # make the x_dagger annotation
    dagger_props = common_arrow_kwargs()
    ax.annotate(xytext=(x_min,y_low_plot),xy=(x_max,y_low_plot),
                s=r"x$^{\ddag}$",**dagger_props)
    # make the delta_G_dagger annotation
    ax.annotate(xytext=(x_max,0),xy=(x_max,y_max),s=r"$\Delta$G$^{\ddag}$",
                **dagger_props)
    # make the extension scale bar 
    width = scale_fraction_width * (max(x)-min(x))
    label = "{:.2g}nm".format(width)
    x_text = x_min
    y_text = np.mean(plt.ylim())
    PlotUtilities.scale_bar_x(x=x_text,
                              y=y_text,
                              s=label,
                              width=width)    
    PlotUtilities.lazyLabel("Extension","Free Energy (k$_B$T)","")
    plt.ylim(-4,max(plt.ylim()))


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    loading_rate_histogram,rupture_forces_histograms,models,rupture_forces,\
        mean_rupture_forces,stdev_rupture_forces,x,landscape = \
            CheckpointUtilities.getCheckpoint("./ruptures.pkl",
                                              generate_rupture_histograms,
                                              True)
    # get the <min.max>(mean + stdev) for the rutpure force (used to set the 
    # intelligently)
    mean_rupture_lower = 0
    mean_rupture_upper = np.max(mean_rupture_forces+stdev_rupture_forces)
    rupture_limits = np.array([mean_rupture_lower,mean_rupture_upper])*1e12
    data_file = "../_Data/example_protein.pkl"
    data = CheckpointUtilities.lazy_load(data_file)
    split_fec,info_final = Detector._predict_full(data)
    # get the plotting versions of the time, etc
    time= split_fec.retract.Time
    time_plot = time - time[0]
    force_plot = split_fec.retract.Force * 1e12
    force_interp_plot = split_fec.retract_spline_interpolator()(time) * 1e12
    # plot everything
    n_cols = 3
    n_rows = 2
    gs = gridspec.GridSpec(n_rows, n_cols)
    ylim_force_pN = [-40,max(force_interp_plot)*1.2]
    ylim_prob = [min(info_final.cdf)/2,2]
    arrow_kwargs = dict(plot_x=time_plot,plot_y=force_plot,
                        markersize=50)
    fig = PlotUtilities.figure(figsize=(16,8))
    # # plot the 'raw' force
    ax1 = plt.subplot(gs[0,0])
    plot_fec_scaled(time_plot,force_plot,force_interp_plot,info_final,
                    arrow_kwargs)
    plt.ylim(ylim_force_pN)
    # # plot the energy landscape with annotations
    ax = plt.subplot(gs[0,1:])
    plot_landscape(x,landscape)
    # # plot the 'zoomed' axis
    ax_zoom = plt.subplot(gs[1,0])
    plot_zoomed(time_plot,force_plot,info_final,ax1,arrow_kwargs)
    # # plot (a single) histogram and model. This one is special, so we 
    # use a slightly different error bar for it 
    fmt_error = dict(marker='v',color='k',markersize=15,linewidth=3)
    example_idx = 4
    loading_rate_example_pN_per_s = loading_rate_histogram[example_idx] * 1e12
    plt.subplot(gs[1,1])
    plot_histogram_and_model(rupture_forces,
                             rupture_forces_histograms[example_idx],
                             models[example_idx],fmt_error)
    PlotUtilities.lazyLabel(r"F$_r$ (pN)","Count","")
    # give the loading rate as an annotation
    plt.text(x=np.mean(rupture_limits),y=np.mean(plt.ylim())*0.5,
             s=df_dt_string+"={:.2g} pN/s".\
             format(loading_rate_example_pN_per_s),fontsize=fontsize,
             horizontalalignment='center',
             verticalalignment='center')
    plt.xlim(rupture_limits)
    plt.subplot(gs[1,2])
    # # plot the distribution of expected rupture forces
    plot_mean_rupture(rupture_forces_histograms,loading_rate_histogram,
                      mean_rupture_forces,stdev_rupture_forces)
    PlotUtilities.lazyLabel(df_dt_string + " (pN/s)",
                            r"<F$_r$> (pN)","")
    # plot the 'extra' point specially, so it sticks out. 
    mean = 1e12*np.mean(rupture_forces_histograms[example_idx])
    stdev = 1e12*np.std(rupture_forces_histograms[example_idx])
    plt.errorbar(x=loading_rate_histogram[example_idx]*1e12,
                 y=mean,
                 yerr=stdev,
                 **fmt_error)
    plt.ylim(rupture_limits)
    PlotUtilities.savefig(fig,"./diagram.png")

if __name__ == "__main__":
    run()
