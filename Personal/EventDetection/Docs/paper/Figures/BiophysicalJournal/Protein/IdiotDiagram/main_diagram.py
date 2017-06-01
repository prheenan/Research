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
from FitUtil.WormLikeChain.Python.Code.WLC_Utils import WlcNonExtensible
from mpl_toolkits.axes_grid.inset_locator import inset_axes


# plotting constants    
raw_force_kwargs = dict(color='k',alpha=0.3)
interp_force_kwargs = dict(color='b',linewidth=1.5)
probabiity_kwargs = dict(color='r')
# how big the scale bars are
scale_fraction_width = 0.13
scale_fraction_offset = 0.3
df_dt_string = r"$\mathbf{\partial}$F/$\mathbf{\partial}$t"
rupture_string = r"F$_{\mathbf{R}}$"
fontsize=8
scale_line_width=1.5

def generate_rupture_histograms():
    num_loading_rates = 10
    np.random.seed(42)
    kbT =  4.1e-21
    delta_G_ddagger = 20 *kbT
    x_ddagger= 0.6e-9
    k0 = 0.1
    nu = 2/3
    beta = 1/kbT
    n_samples = 1000
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
    x = np.linspace(start=-x_ddagger*1.0,stop=x_ddagger*0.8)
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
    style_interp_after['color'] = 'r'
    style_interp_final = dict(**interp_force_kwargs)
    style_interp_final['color'] = 'g'
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
                                      fudge_y=30,**arrow_kwargs)
    # add a scale bar
    PlotUtilities.no_x_label(plt.gca())               
    max_time = max(time_plot)
    width = scale_fraction_width * max_time
    scale_bar_dict = dict(y_label_frac=0.2,linewidth=scale_line_width)
    PlotUtilities.x_scale_bar_and_ticks(scale_bar_dict=scale_bar_dict)

def common_text_kwargs():
    return dict(xycoords='data',
                textcoords='data',
                verticalalignment='center',
                horizontalalignment='center')

def common_arrow_kwargs(arrowprops=dict(arrowstyle="<->",shrinkA=0,
                                        shrinkB=0,linewidth=1,
                                        connectionstyle="arc3")):
    return  dict(fontsize=fontsize,
                 arrowprops=arrowprops,**common_text_kwargs())

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
    _,_,loading_rate,rupture_force,_ = \
            Analysis._loading_rate_helper(x_fit,y_fit)
    raw_kwargs = dict(**raw_force_kwargs)
    raw_kwargs['alpha'] = 0.5
    plt.plot(time_slice,force_slice,**raw_kwargs)
    color_loading_line = 'b'
    xlim = [time_slice[0],time_slice[-1]]
    plt.xlim(xlim)
    ax_zoom = plt.gca()
    PlotUtilities.zoom_effect01(ax1, ax_zoom, *xlim,linestyle=':',
                                color='m')
    # plot annotations showing the loading rate
    bbox_props = dict(boxstyle="rarrow,pad=0.3",linestyle='--',
                      fc=color_loading_line,ec=color_loading_line,
                      alpha=0.3,lw=2)
    fudge = (max(x_fit)-min(x_fit)) * 0.2
    min_time = min(x_fit) + fudge
    event_force = rupture_force
    event_time = time_plot[event_zoom]
    arrow_kw_thick = common_arrow_kwargs()
    linewidth_common = 1.5
    arrow_kw_thick['arrowprops']['linewidth'] = linewidth_common
    arrow_kw_thick['arrowprops']['arrowstyle'] = '->'
    zoom_event_only = [event_zoom]
    # get the x locations of two predicted events
    predicted_x_fracs = [0.8,0.95]
    xlim = plt.xlim()
    ylim = plt.ylim()
    predicted_x = [xlim[0] + (xlim[1]-xlim[0]) * r for r in predicted_x_fracs]
    # plot them as lines
    style_predicted = dict(linestyle=':',color='c',linewidth=linewidth_common)
    style_arrow_predicted = dict(**arrow_kw_thick)
    style_arrow_predicted['color'] = 'b'
    for x in predicted_x:
        plt.axvline(x,**style_predicted)
    text_box_kwargs = dict(horizontalalignment='center',fontweight='bold',
                           verticalalignment='center',fontsize=fontsize)
    # plot a text box on top of the lines
    predicted_line_style = ":"
    plt.text(x=np.mean(predicted_x),
             y=np.mean(ylim)+abs(np.diff(ylim))*0.1,backgroundcolor='c',
             bbox=dict(linestyle=predicted_line_style,color='c'),color='w',
             s="Predicted events",rotation=90,zorder=10,**text_box_kwargs)
    # plot a text box for the event
    plt.axvline(event_time,color='g')
    plt.text(x=event_time-abs(np.diff(plt.xlim()))*0.05,
             y=np.mean(ylim),backgroundcolor='g',
             bbox=dict(linestyle='-',color='g'),color='w',
             s="True event",rotation=90,zorder=10,**text_box_kwargs)
    # plot an from the event to the closest x 
    closest_x_idx = np.argmin(np.abs(predicted_x-event_time))
    closest_x = predicted_x[closest_x_idx]
    dx = closest_x-event_time
    dy = 0
    plt.text(event_time+dx/2,event_force*0.95,r"d$_{t\rightarrow p}$",
             color='g',**text_box_kwargs)
    plt.annotate("", xytext=(event_time, event_force), 
                 xy=(closest_x, event_force),
                 arrowprops=dict(arrowstyle="->",color='g',
                                 linewidth=linewidth_common))
    # plot arrows from the predictions to the actual
    # see:
    # matplotlib.org/api/pyplot_api.html?highlight=arrow#matplotlib.pyplot.arrow
    # we need to use arowstyle 
    ax = plt.gca()
    y_text= np.mean(ylim)*1.07
    plot_y = [y_text*1.07,y_text*1.12]
    for x,plot_y_tmp in zip(predicted_x,plot_y):
        ax.annotate("", xy=(event_time, plot_y_tmp), xytext=(x, plot_y_tmp),
                    arrowprops=dict(arrowstyle="->",color='c',
                                    linewidth=linewidth_common,
                                    linestyle='-'))
    plt.text(event_time+dx/2,y_text,r"d$_{p\rightarrow t}$",
             color='c',**text_box_kwargs)
    # add a scalebar...
    scale_bar_kwargs = dict(dict(y_frac=0.5,y_label_frac=0.1,
                                 linewidth=scale_line_width))
    PlotUtilities.x_scale_bar_and_ticks(scale_bar_kwargs)


def plot_mean_rupture(rupture_forces_histograms,loading_rate_histogram,
                      mean_rupture_forces,stdev_rupture_forces):
    sizes = [r.size for r in rupture_forces_histograms]
    loading_plot = loading_rate_histogram*1e12
    mean_plot = mean_rupture_forces*1e12                                
    stdev_plot = stdev_rupture_forces*1e12
    sem_plot = stdev_plot/sizes
    # plot the theory first
    mean_rupture_string = "<" + rupture_string + ">"
    plt.plot(loading_plot,mean_plot,linewidth=1.5)
    # plot the data on top of the theory
    means = np.mean(rupture_forces_histograms,axis=1) *1e12
    stdevs= np.std(rupture_forces_histograms,axis=1) * 1e12
    sems = stdevs/sizes 
    plt.plot(loading_plot,means,'ro',
             markersize=2)
    plt.gca().set_xscale('log')


def plot_histogram_and_model(rupture_forces,rupture,model,kwargs_errorbar,
                             loading_rate_pN_per_s,bins=25):
    rupture_plot = rupture*1e12
    n,_,_ = plt.hist(rupture_plot,alpha=0.5,bins=bins,edgecolor='k',linewidth=1)
    plt.plot(rupture_forces*1e12,model*max(n)/max(model),**kwargs_errorbar)
    plt.xlim([min(rupture_plot),max(rupture_plot)])
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.ylim([ylim[0],ylim[1]*1.2])



def plot_landscape(x,landscape,ax=plt.gca()):
    landscape -= min(landscape)
    landscape /= 4.1e-21
    x *= 1e9
    plt.plot(x,landscape,color='k',linewidth=1.5)
    PlotUtilities.no_x_label(ax)
    # determine where to put all the annotations
    min_idx =np.argmin(landscape)
    max_idx = min_idx + np.argmax(landscape[min_idx:])
    x_max = x[max_idx]
    x_min = x[min_idx]
    x_range = x_max-x_min
    y_low_plot = -3
    y_max = landscape[max_idx]
    y_range = y_max - landscape[min_idx]
    # make the x_dagger annotation
    dagger_props = common_arrow_kwargs()
    text_kwargs = dict(horizontalalignment='center',
                       verticalalignment='center',fontsize=8)
    dagger_props['arrowprops']['arrowstyle'] = '<->'
    ax.annotate(xytext=(x_min,y_low_plot),xy=(x_max,y_low_plot),
                s=r"",**dagger_props)
    ax.text(x=np.mean([x_max,x_min]),y=y_low_plot+y_range*0.1,s="x$^{\ddag}$",
            **text_kwargs)
    # make the delta_G_dagger annotation
    ax.annotate(xytext=(x_max,0),xy=(x_max,y_max),s=r"",**dagger_props)
    ax.text(x=x_max+x_range*0.4,y=np.mean(landscape)*0.7,
            s="$\Delta$G$^{\ddag}$",**text_kwargs)
    # make the k0 annotation
    fudge = (x_max-x_min)*0.5
    y_k0 = y_max * 0.95
    ax.annotate('',
                xy=(x_max+fudge, y_k0), xycoords='data',
                xytext=(x_max-fudge, y_k0), textcoords='data',
                arrowprops=dict(arrowstyle="->",shrinkA=0,shrinkB=0,
                                connectionstyle="angle3,angleA=45,angleB=-45"))
    ax.text(x=x_max,y=y_k0*1.25,s="k$_0$",**text_kwargs)
    plt.ylim(y_low_plot*2,max(plt.ylim())*1.1)
    xlim = plt.xlim()
    plt.xlim(xlim[0],xlim[1]+x_range * 0.4)
    PlotUtilities.no_y_ticks(ax=ax)
    PlotUtilities.no_x_ticks(ax=ax)
    PlotUtilities.no_y_label(ax=ax)
    PlotUtilities.no_x_label(ax=ax)

def cantilever_image_plot(image_location):
    in_ax = plt.gca()
    im = plt.imread(image_location)
    in_ax.imshow(im,interpolation="bicubic",aspect='auto')
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
    lazy_kwargs = dict()
    data_file = "../_Data/example_protein.pkl"
    image_location =  "../_Data/JILA_PFC_PRH_Edi.png"
    data = CheckpointUtilities.lazy_load(data_file)
    split_fec,info_final = Detector._predict_full(data)
    # get the plotting versions of the time, etc
    time= split_fec.retract.Time
    time_plot = time - time[0]
    force_plot = split_fec.retract.Force * 1e12
    force_interp_plot = split_fec.retract_spline_interpolator()(time) * 1e12
    # plot everything
    n_cols = 2
    n_rows = 3
    # 'master' grid spec is 2x1
    gs0 = gridspec.GridSpec(2,1)
    gs = gridspec.GridSpecFromSubplotSpec(n_rows, n_cols, subplot_spec=gs0[0],
                                          hspace=0.75,wspace=0.5)
    ylim_force_pN = [-35,max(force_interp_plot)*1.3]
    ylim_prob = [min(info_final.cdf)/2,2]
    arrow_kwargs = dict(plot_x=time_plot,plot_y=force_plot,
                        markersize=75)
    fig = PlotUtilities.figure(figsize=(3.25,6))
    # # plot the experimental image
    in_ax = plt.subplot(gs[:,0])
    cantilever_image_plot(image_location)
    # remove its border
    for spine in in_ax.spines.values():
        spine.set_visible(False)
    # # plot the cartoon of the fec; easier to just call out the axis
    ax_fec = plt.axes([0.635,0.86,0.31,0.105])
    # define the regions where we are attached to the molecule...
    slice_tuples = [ [slice(0,0.05),0.08],
                     [slice(0,0.2),0.25],
                     [slice(0,0.4),0.5],
                     [slice(0,0.6),0.7],
                     [slice(0,0.8),0.92],
                     [slice(0,0.95),1.05]]
    n_cartoon = 500
    x_cartoon = np.linspace(0,1.0,num=n_cartoon)
    f = np.zeros(x_cartoon.size)
    slice_starts = [0]
    for slice_active,L0 in slice_tuples:
        slice_tmp = slice(int(slice_active.start * n_cartoon),
                          int(slice_active.stop * n_cartoon))
        slice_starts.append(slice_tmp.stop)
        f[slice_tmp] += WlcNonExtensible(ext=x_cartoon[slice_tmp],
                                         kbT=1,Lp=max(x)/30,L0=L0)
    colors = ['b','r','r','r','r','g']
    slice_starts.append(None)
    for i in range(len(slice_starts)-2):
        slice_tmp = slice(slice_starts[i],slice_starts[i+1]+1,1)
        plt.plot(x_cartoon[slice_tmp],f[slice_tmp],color=colors[i])
    PlotUtilities.lazyLabel("Time","Force","Analysis scheme",**lazy_kwargs )
    PlotUtilities.no_x_label(ax_fec)
    PlotUtilities.no_x_ticks(ax_fec)
    PlotUtilities.no_y_label(ax_fec)
    PlotUtilities.no_y_ticks(ax_fec)
    # # plot the loading rate stuff
    fmt_error = dict(linewidth=1.5)
    example_idx = 4
    loading_rate_example_pN_per_s = loading_rate_histogram[example_idx] * 1e12
    plt.subplot(gs[1,1])
    plot_histogram_and_model(rupture_forces,
                             rupture_forces_histograms[example_idx],
                             models[example_idx],fmt_error,
                             loading_rate_example_pN_per_s)
    PlotUtilities.lazyLabel(rupture_string+ " (pN)","Count","",**lazy_kwargs)
    PlotUtilities.tom_ticks(num_major=3,change_x=False)
    plt.subplot(gs[2,1])
    # # plot the distribution of expected rupture forces
    plot_mean_rupture(rupture_forces_histograms,loading_rate_histogram,
                      mean_rupture_forces,stdev_rupture_forces)
    PlotUtilities.lazyLabel(df_dt_string + " (pN/s)",
                            rupture_string + " (pN)","",**lazy_kwargs)
    plt.ylim(rupture_limits)
    # # plot the energy landscape with annotations
    # add axes is [left,bottom,width,height]
    in_ax_landscape= fig.add_axes([0.085,0.72,0.225,0.17])
    plot_landscape(x,landscape,ax=in_ax_landscape)
    energy_kwargs = dict(axis_kwargs=dict(fontsize=8))
    # remove the upper and right part of the frames
    in_ax_landscape.spines['right'].set_visible(False)
    in_ax_landscape.spines['top'].set_visible(False)
    PlotUtilities.lazyLabel("Extension","Free Energy","",**energy_kwargs)
    # # # Second gridspec (se we can more easily control wasted space)
    gs_data = gridspec.GridSpecFromSubplotSpec(3,1, 
                                               subplot_spec=gs0[1])
    # # plot the 'raw' force
    ax1 = plt.subplot(gs_data[0,:])
    plot_fec_scaled(time_plot,force_plot,force_interp_plot,info_final,
                    arrow_kwargs)
    PlotUtilities.lazyLabel("","F (pN)",
                            "Extracting rupture properties",
                            **lazy_kwargs)
    PlotUtilities.tom_ticks(num_major=5,change_x=False)
    xlim = plt.xlim()
    plt.ylim(ylim_force_pN)
    # # plot the 'zoomed' axis
    ax_zoom = plt.subplot(gs_data[1:,:])
    plot_zoomed(time_plot,force_plot,info_final,ax1,arrow_kwargs)
    PlotUtilities.lazyLabel("Time","F (pN)","",**lazy_kwargs)
    PlotUtilities.no_x_label(ax_zoom)
    PlotUtilities.tom_ticks(num_major=5,change_x=False)
    axis_func = lambda axes: [a for i,a in enumerate(axes) if i != 4]
    loc_subplot = [-0.5,1.17]
    locs = [ [-0.45,1.02],
             loc_subplot,
             loc_subplot,
             loc_subplot,
             [-0.18,1.15],
             [-0.18,0.95]]
    PlotUtilities.label_tom(fig,axis_func=axis_func,loc=locs)
    PlotUtilities.savefig(fig,"./diagram.png")

if __name__ == "__main__":
    run()
