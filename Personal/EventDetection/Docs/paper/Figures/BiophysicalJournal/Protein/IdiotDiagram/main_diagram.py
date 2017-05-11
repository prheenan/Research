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
from mpl_toolkits.axes_grid.inset_locator import inset_axes


# plotting constants    
raw_force_kwargs = dict(color='k',alpha=0.3)
interp_force_kwargs = dict(color='b',linewidth=3)
probabiity_kwargs = dict(color='r')
# how big the scale bars are
scale_fraction_width = 0.13
scale_fraction_offset = 0.3
df_dt_string = r"dF/dt"
rupture_string = r"F$_{\mathrm{r}}$"
fontsize=18

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
    # get the style for before and after
    style_interp_before = dict(**interp_force_kwargs)
    style_interp_before['color'] = 'b'
    style_interp_after = dict(**interp_force_kwargs)
    style_interp_after['color'] = 'g'
    style_interp_final = dict(**interp_force_kwargs)
    style_interp_final['color'] = 'm'
    # get the slices
    event_final = info_final.event_idx[-2]
    event_initial = info_final.event_idx[0]
    slice_before = slice(0,event_initial,1)
    slice_after = slice(event_initial,event_final,1)
    slice_final = slice(event_final,None,1)
    plt.plot(time_plot,force_plot,**raw_force_kwargs)
    plt.plot(time_plot[slice_before],force_interp_plot[slice_before],
             **style_interp_before)
    plt.plot(time_plot[slice_after],force_interp_plot[slice_after],
             **style_interp_after)
    plt.plot(time_plot[slice_final],force_interp_plot[slice_final],
             **style_interp_final)
    PlotUtilities.lazyLabel("","Force (pN)","")
    # plot arrows above the events
    Plotting.plot_arrows_above_events(event_idx=info_final.event_idx,
                                      fudge_y=20,**arrow_kwargs)
    # add a scale bar
    PlotUtilities.no_x_label(plt.gca())               
    max_time = max(time_plot)
    width = scale_fraction_width * max_time
    label = "{:.1g}s".format(width)
    PlotUtilities.scale_bar_x(x=scale_fraction_offset*max_time,
                              y=np.max(plt.ylim())*0.7,s=label,
                              width=width,fontsize=fontsize)  

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
    n_around_idx = 2e-3
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
    rotation =20
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
    arrow_kw_thick['arrowprops']['arrowstyle'] = '->'
    ax_zoom.annotate(rupture_string,
                     xy=(event_time,event_force), 
                     xytext=(min_time, event_force),**arrow_kw_thick)
    zoom_event_only = [event_zoom]
    # make the marker larger for this subplot
    arrow_kwargs = dict(**arrow_kwargs)
    arrow_kwargs['markersize'] *= 3
    Plotting.plot_arrows_above_events(event_idx=zoom_event_only,fudge_y=6,
                                      **arrow_kwargs)
    # get the x locations of two predicted events
    predicted_x_fracs = [0.8,0.95]
    xlim = plt.xlim()
    ylim = plt.ylim()
    predicted_x = [xlim[0] + (xlim[1]-xlim[0]) * r for r in predicted_x_fracs]
    # plot them as lines
    style_predicted = dict(linestyle='-.',color='c',linewidth=3)
    style_arrow_predicted = dict(**arrow_kw_thick)
    style_arrow_predicted['color'] = 'b'
    for x in predicted_x:
        plt.axvline(x,**style_predicted)
    text_box_kwargs = dict(horizontalalignment='center',
                           verticalalignment='center',fontsize=fontsize)
    # plot a text box on top of the lines
    plt.text(x=np.min(predicted_x)+abs(np.diff(predicted_x))*0.2,
             y=np.mean(ylim),backgroundcolor='c',
             bbox=dict(linestyle='-.',color='c'),color='w',
             s="Predictions",rotation=90,zorder=10,**text_box_kwargs)
    # plot an from the event to the closest x 
    closest_x_idx = np.argmin(np.abs(predicted_x-event_time))
    closest_x = predicted_x[closest_x_idx]
    dx = closest_x-event_time
    dy = 0
    plt.text(event_time+dx/2,event_force*1.05,r"d$_{t\rightarrow p}$",
             color='g',**text_box_kwargs)
    plt.annotate("", xytext=(event_time, event_force), 
                 xy=(closest_x, event_force),
                arrowprops=dict(arrowstyle="->",color='g',linewidth=2))
    # plot arrows from the predictions to the actual
    # see:
    # matplotlib.org/api/pyplot_api.html?highlight=arrow#matplotlib.pyplot.arrow
    # we need to use arowstyle 
    ax = plt.gca()
    y_text= np.mean(ylim)*1.07
    plot_y = [y_text*1.07,y_text*1.12]
    for x,plot_y_tmp in zip(predicted_x,plot_y):
        ax.annotate("", xy=(event_time, plot_y_tmp), xytext=(x, plot_y_tmp),
                    arrowprops=dict(arrowstyle="->",color='c',linewidth=3,
                                    linestyle="-."))
    plt.text(event_time+dx/2,y_text,r"d$_{p\rightarrow t}$",
             color='c',**text_box_kwargs)
    # add a scalebar...
    dx_zoom_full =abs(time_slice[-1]-time_slice[0])
    width = scale_fraction_width * dx_zoom_full
    label = "{:.1g}ms".format(1000*width)
    x_text = time_slice[0] + dx_zoom_full*scale_fraction_offset
    PlotUtilities.scale_bar_x(x=x_text,
                              y=min(force_slice)*1.1,s=label,
                              width=width,fontsize=fontsize)    
    PlotUtilities.lazyLabel("Time","Force (pN)","") 
    PlotUtilities.no_x_label(ax)


def plot_mean_rupture(rupture_forces_histograms,loading_rate_histogram,
                      mean_rupture_forces,stdev_rupture_forces):
    sizes = [r.size for r in rupture_forces_histograms]
    loading_plot = loading_rate_histogram*1e12
    mean_plot = mean_rupture_forces*1e12                                
    stdev_plot = stdev_rupture_forces*1e12
    sem_plot = stdev_plot/sizes
    # plot the theory first
    mean_rupture_string = "<" + rupture_string + ">"
    plt.plot(loading_plot,mean_plot,
             label=(mean_rupture_string + " (Theory)"),
             linewidth=3)
    stdev_style = dict(linestyle=":",linewidth=3,color='m')
    plt.plot(loading_plot,mean_plot-sem_plot,
             **stdev_style)
    plt.plot(loading_plot,mean_plot+sem_plot,**stdev_style)
    # plot the data on top of the theory
    means = np.mean(rupture_forces_histograms,axis=1) *1e12
    stdevs= np.std(rupture_forces_histograms,axis=1) * 1e12
    sems = stdevs/sizes 
    plt.errorbar(loading_plot,y=means,yerr=sems,fmt='ro')
    plt.gca().set_xscale('log')


def plot_histogram_and_model(rupture_forces,rupture,model,kwargs_errorbar):
    rupture_plot = rupture*1e12
    n,_,_ = plt.hist(rupture_plot,alpha=0.5)
    plt.plot(rupture_forces*1e12,model*max(n)/max(model),linewidth=3)

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
    dagger_props['arrowprops']['arrowstyle'] = '<->'
    ax.annotate(xytext=(x_min,y_low_plot),xy=(x_max,y_low_plot),
                s=r"x$^{\ddag}$",**dagger_props)
    # make the delta_G_dagger annotation
    ax.annotate(xytext=(x_max,0),xy=(x_max,y_max),s=r"$\Delta$G$^{\ddag}$",
                **dagger_props)
    # make the extension scale bar 
    width = 0.2 * (max(x)-min(x))
    label = "{:.1g}nm".format(width)
    x_text = x_min
    y_text = np.mean(plt.ylim())
    PlotUtilities.scale_bar_x(x=x_text,
                              y=y_text,
                              s=label,fontsize=fontsize,
                              width=width)    
    PlotUtilities.lazyLabel("Extension","Free Energy (k$_B$T)","")
    plt.ylim(-4,max(plt.ylim()))

def cantilever_image_plot(image_location):
    in_ax = plt.gca()
    im = plt.imread(image_location)
    in_ax.imshow(im,interpolation="bicubic")
    in_ax.set_xticks([])
    in_ax.set_yticks([])


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
    image_location =  "../_Data/pulling_figure_4nug2.png"
    data = CheckpointUtilities.lazy_load(data_file)
    split_fec,info_final = Detector._predict_full(data)
    # get the plotting versions of the time, etc
    time= split_fec.retract.Time
    time_plot = time - time[0]
    force_plot = split_fec.retract.Force * 1e12
    force_interp_plot = split_fec.retract_spline_interpolator()(time) * 1e12
    # plot everything
    n_cols = 3
    n_rows = 5
    gs = gridspec.GridSpec(n_rows, n_cols)
    ylim_force_pN = [-35,max(force_interp_plot)*1.2]
    ylim_prob = [min(info_final.cdf)/2,2]
    arrow_kwargs = dict(plot_x=time_plot,plot_y=force_plot,
                        markersize=75,marker=u'$\u2193$')
    fig = PlotUtilities.figure(figsize=(8,11))
    # # plot the experimental image
    in_ax = plt.subplot(gs[:3,0])
    cantilever_image_plot(image_location)
    # # plot the 'raw' force
    ax1 = plt.subplot(gs[3,:])
    plot_fec_scaled(time_plot,force_plot,force_interp_plot,info_final,
                    arrow_kwargs)
    xlim = plt.xlim()
    plt.ylim(ylim_force_pN)
    # # plot the 'zoomed' axis
    ax_zoom = plt.subplot(gs[4,:])
    plot_zoomed(time_plot,force_plot,info_final,ax1,arrow_kwargs)
    # # plot (a single) histogram and model. This one is special, so we 
    # use a slightly different error bar for it 
    fmt_error = dict(marker='v',color='k',markersize=10,linewidth=3)
    example_idx = 4
    loading_rate_example_pN_per_s = loading_rate_histogram[example_idx] * 1e12
    plt.subplot(gs[1,1:])
    plot_histogram_and_model(rupture_forces,
                             rupture_forces_histograms[example_idx],
                             models[example_idx],fmt_error)
    PlotUtilities.lazyLabel(rupture_string+ " (pN)","Count","")
    # give the loading rate as an annotation
    loading_rate_str = ("Simulation\n" + r"$\frac{dF}{dt}$" + "={:.2g} pN/s".\
                        format(loading_rate_example_pN_per_s))
    plt.text(x=np.mean(rupture_limits),y=np.mean(plt.ylim())*0.6,
             s=loading_rate_str,fontsize=fontsize,
             horizontalalignment='center',
             verticalalignment='center')
    plt.xlim(rupture_limits)
    plt.subplot(gs[2,1:])
    # # plot the distribution of expected rupture forces
    plot_mean_rupture(rupture_forces_histograms,loading_rate_histogram,
                      mean_rupture_forces,stdev_rupture_forces)
    PlotUtilities.lazyLabel(df_dt_string + " (pN/s)",
                            rupture_string + " (pN)","",
                            legend_kwargs=dict(handlelength=1))
    plt.ylim(rupture_limits)
    # # plot the energy landscape with annotations
    ax = plt.subplot(gs[0,1:])
    plot_landscape(x,landscape)
    PlotUtilities.savefig(fig,"./diagram.png",
                          subplots_adjust=dict(hspace=0.4))

if __name__ == "__main__":
    run()
