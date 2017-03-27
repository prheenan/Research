# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,string

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection._2SplineEventDetector import Detector

from Research.Personal.EventDetection.Util import Analysis,Plotting

import matplotlib.gridspec as gridspec

# style keywords for this file only
common = dict(rasterized=True)
style_approach = dict(color='k', **common)
style_raw = dict(alpha=0.3,**common)
style_filtered = dict(linewidth=3,**common)
retract_style = dict(color='g',**common)
style_retract_error_dist = dict(linewidth=2,**common)
label_s_t = r"s$_t$"
label_r_t = r"r$_t$"
lazy_kwargs = dict(loc='lower right',frameon=False)

error_label = "r$_\mathrm{t}$ or\ns$_\mathrm{t}$ (pN)"
filtered_error_label = "s$_{\mathrm{t}}$$^{*}$ (pN)"
force_label = "Force (pN)"

def tick_function():
    tick_kwargs = dict(num_x_major=4,num_y_major=4)
    PlotUtilities.tick_axis_number(**tick_kwargs)


def plot_fec(time_approach,force_approach,interp_approach,
             ylim_force,xlim_approach,label):
    plt.plot(time_approach,force_approach,label=label,
             alpha=0.3,**style_approach)
    plt.plot(time_approach,interp_approach,label="Spline ($g^{*}_t$)))",
             **style_approach)
    PlotUtilities.lazyLabel("",force_label,
                            "Estimating $\epsilon$ and $\sigma$",
                            frameon=False)
    plt.xlim(xlim_approach)
    plt.ylim(ylim_force)

def plot_retract_fec(x_plot,force_plot,slice_before,slice_after,
                     force_filtered_plot,ylim_force):
    Plotting.before_and_after(x_plot,force_plot,slice_before,slice_after,
                              style_raw,label="Raw (retract)")
    Plotting.before_and_after(x_plot,force_filtered_plot,
                              slice_before,slice_after,style_filtered,
                              label="Spline ($g^{*}_t$)))")
    title = "Calculating the no-event probability"
    plt.ylim(ylim_force)
    PlotUtilities.lazyLabel("","",title,**lazy_kwargs)

def plot_approach_error(time_approach_plot,approach_diff,approach_stdev_plot,
                        xlim_approach,ylim_diff):
    plt.plot(time_approach_plot,approach_diff,alpha=0.3,
             label=label_r_t,**style_approach)
    plt.plot(time_approach_plot,approach_stdev_plot,label=label_s_t,
             linestyle='--',**style_approach)
    plt.xlim(xlim_approach)
    PlotUtilities.no_x_label()
    tick_function()
    plt.ylim(*ylim_diff)
    PlotUtilities.lazyLabel("",error_label,"",frameon=True,loc='upper right')

def plot_retract_error(x_plot,diff_pN,slice_before,slice_after,stdev_plot,
                       ylim_diff):
    filted_stdev_style = dict(linestyle='--',**style_filtered)
    Plotting.before_and_after(x_plot,diff_pN,slice_before,slice_after,
                              style_raw,label=label_r_t)
    Plotting.before_and_after(x_plot,stdev_plot,slice_before,slice_after,
                              filted_stdev_style,label=label_s_t)
    PlotUtilities.no_x_label()
    plt.ylim(*ylim_diff)
    tick_function()
    PlotUtilities.lazyLabel("","","",frameon=True,loc='upper right')


def plot_filtered_stdev(time_approach_plot,approach_stdev_plot,
                        epsilon_plot,sigma_plot,xlim_approach,
                        ylim_diff_filtered):
    plt.plot(time_approach_plot,approach_stdev_plot,**style_approach)
    plot_epsilon(epsilon_plot,sigma_plot)
    plt.xlim(xlim_approach)
    plt.ylim(*ylim_diff_filtered)
    tick_function()
    PlotUtilities.lazyLabel("Time (s)",filtered_error_label,"",
                            frameon=False,loc='upper right')


def plot_filtered_retract_stdev(x_plot,stdev_plot,slice_before,slice_after,
                                epsilon_plot,sigma_plot,ylim_diff_filtered):
    Plotting.before_and_after(x_plot,stdev_plot,slice_before,slice_after,
                              style_retract_error_dist)
    plot_epsilon(epsilon_plot,sigma_plot)
    PlotUtilities.lazyLabel("","","",frameon=False,
                            loc="upper right")
    PlotUtilities.no_x_label()
    plt.ylim(*ylim_diff_filtered)
    tick_function()

def plot_probability(threshold,x_plot,prob_final,slice_before,slice_after):
    plt.yscale('log')
    plt.axhline(threshold,linewidth=3,linestyle='--',
                label="threshold",color='k')
    Plotting.before_and_after(x_plot,prob_final,slice_before,slice_after,
                              style_retract_error_dist)
    PlotUtilities.lazyLabel("Time (s)","Probability","",**lazy_kwargs)
    plt.ylim([min(plt.ylim()),1.5])



def plot_epsilon(epsilon_plot,sigma_plot):
    epsilon_style = dict(color='b')
    plt.axhline(epsilon_plot,label="$\epsilon$")
    plt.axhline(epsilon_plot+sigma_plot,linestyle='--',
                label="$\epsilon \pm \sigma$",**epsilon_style)
    plt.axhline(epsilon_plot-sigma_plot,linestyle='--',**epsilon_style)

def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    out_fig = "algorithm.pdf"
    force = False
    example = read_and_cache_file(data_base + "rupture.csv",has_events=True,
                                  force=force,cache_directory=data_base)
    n_filter = 1000
    kw = dict(cache_directory=data_base,force=force)
    events = example.Events
    fec_split = Analysis.zero_and_split_force_extension_curve(example)
    event_idx_retract = fec_split.get_retract_event_centers()
    event_idx = event_idx_retract[0]
    zoom_factor = 20
    points_around = int(np.ceil(event_idx/zoom_factor))
    retract = fec_split.retract
    retract.Force -= np.median(retract.Force)
    # get everything in terms of ploting variables
    x_plot_f = lambda x: x - min(x)
    y_plot_f = lambda y: y * 1e12
    time = retract.Time
    force = retract.Force
    # get the slices for the epsilon and sigmas
    n = fec_split.approach.Force.size
    offset_approach =  n-fec_split.get_predicted_approach_surface_index()
    slice_fit_approach = slice(offset_approach,-offset_approach,1)
    # get the approach 
    approach = fec_split.approach
    approach_time = approach.Time
    time_approach = x_plot_f(approach_time)
    force_approach = y_plot_f(approach.Force)
    interpolator_approach = fec_split.approach_spline_interpolator()
    interp_approach = y_plot_f(interpolator_approach(approach_time))
    approach_diff = force_approach-interp_approach
    approach_stdev,_,_ =  fec_split.\
    _approach_stdevs_epsilon_and_sigma(slice_fit_approach=slice_fit_approach)
    # get the retract 
    interpolator = fec_split.retract_spline_interpolator()
    force_filtered = interpolator(time)
    n_points = fec_split.tau_num_points
    diff = force-force_filtered
    diff_pN = y_plot_f(diff)
    stdev,_,_ = Analysis.stdevs_epsilon_sigma(force,force_filtered,
                                              n_points)
    epsilon,sigma = fec_split.calculate_epsilon_and_sigma(n_points=n_points,\
                                    slice_fit_approach=slice_fit_approach)
    # get the probability
    prob,_= Detector._no_event_probability(time,interp=interpolator,
                                           y=force,
                                           n_points=n_points,epsilon=epsilon,
                                           sigma=sigma)
    # get the prediction info
    threshold = 0.1
    _,predict_info = Detector._predict_full(example,threshold=threshold)
    # get the final masks
    mask_final = predict_info.condition_results[-1]
    bool_final = np.zeros(time.size)
    bool_final[mask_final] = 1
    prob_final = predict_info.probabilities[-1]
    # plot everything
    x_plot = time
    xlim_approach = [min(approach_time),min(time)]
    xlim_retract = [min(time),max(time)]
    force_plot = y_plot_f(force)
    force_filtered_plot = y_plot_f(force_filtered)
    epsilon_plot,sigma_plot = y_plot_f(epsilon),y_plot_f(sigma)
    stdev_plot = y_plot_f(stdev)
    before,after = Analysis.get_before_and_after_and_zoom_of_slice(fec_split)
    slice_before,slice_after = before[0],after[0]
    time_approach_plot = time_approach[slice_fit_approach]
    approach_stdev_plot = y_plot_f(approach_stdev)
    ylim_diff_filtered = [min(stdev_plot),max(stdev_plot)]
    ylim_diff = [min(diff_pN),max(diff_pN)]
    ylim_force_min = min([min(force_plot),min(force_approach)])
    ylim_force_max = max([max(force_plot),max(force_approach)])
    ylim_force = [ylim_force_min,ylim_force_max]
    fig = PlotUtilities.figure((16,12))
    n_rows = 3
    n_cols = 2
    """
    make just the retract figures
    """
    fig = PlotUtilities.figure((8,10))
    plt.subplot(4,1,1)
    plot_retract_fec(x_plot,force_plot,slice_before,slice_after,
                     force_filtered_plot,ylim_force)
    # remove the title 
    plt.title("")
    PlotUtilities.xlabel("Time (s)")
    PlotUtilities.ylabel(force_label)
    PlotUtilities.legend(frameon=True)
    PlotUtilities.x_label_on_top()
    plt.subplot(4,1,2)
    plot_retract_error(x_plot,diff_pN,slice_before,slice_after,stdev_plot,
                       ylim_diff)
    PlotUtilities.ylabel(error_label)
    plt.subplot(4,1,3)
    plot_filtered_retract_stdev(x_plot,stdev_plot,slice_before,slice_after,
                                epsilon_plot,sigma_plot,ylim_diff_filtered)
    PlotUtilities.ylabel(filtered_error_label)
    PlotUtilities.legend(frameon=True)
    plt.subplot(4,1,4)
    plot_probability(threshold,x_plot,prob_final,slice_before,slice_after)
    PlotUtilities.no_x_label()
    PlotUtilities.xlabel("")
    PlotUtilities.savefig(fig,out_fig.replace(".pdf","_retract.pdf"))
    """
    make just the approach figure
    """
    fig = PlotUtilities.figure((6,7))
    plt.subplot(3,1,1)
    plot_fec(time_approach,force_approach,interp_approach,
             ylim_force,xlim_approach,"Raw (approach)")
    # remove the title 
    PlotUtilities.xlabel("Time (s)")
    PlotUtilities.x_label_on_top()
    plt.title("")
    plt.subplot(3,1,2)
    plot_approach_error(time_approach_plot,approach_diff[slice_fit_approach],
                        approach_stdev_plot,xlim_approach,ylim_diff)
    plt.subplot(3,1,3)
    plot_filtered_stdev(time_approach_plot,approach_stdev_plot,
                        epsilon_plot,sigma_plot,xlim_approach,
                        ylim_diff_filtered)
    tick_function()
    PlotUtilities.no_x_label()
    PlotUtilities.xlabel("")
    # increase the y limits to more or less the distribution
    plt.ylim([ epsilon_plot-sigma_plot*10,epsilon_plot+sigma_plot*10])
    PlotUtilities.savefig(fig,out_fig.replace(".pdf","_approach.pdf"))
    """
    make the plot for the paper
    """
    fig = PlotUtilities.figure((16,8))
    gs = gridspec.GridSpec(4, 2)
    plt.subplot(gs[0,0])
    plot_fec(time_approach,force_approach,interp_approach,
             ylim_force,xlim_approach,"Raw (approach)")
    tick_function()
    PlotUtilities.no_x_label()
    plt.subplot(gs[0,1])
    # plot the 'raw' error distribution for the approach
    plot_retract_fec(x_plot,force_plot,slice_before,slice_after,
                     force_filtered_plot,ylim_force)
    tick_function()
    PlotUtilities.no_x_label()
    plt.subplot(gs[1,0])
    # plot the 'raw error distribution for the approach
    plot_approach_error(time_approach_plot,approach_diff[slice_fit_approach],
                        approach_stdev_plot,xlim_approach,ylim_diff)
    plt.subplot(gs[1,1])
    # filtered error distribution for the retract
    plot_retract_error(x_plot,diff_pN,slice_before,slice_after,stdev_plot,
                       ylim_diff)
    plt.subplot(gs[2,0])
    plot_filtered_stdev(time_approach_plot,approach_stdev_plot,
                        epsilon_plot,sigma_plot,xlim_approach,
                        ylim_diff_filtered)
    # filtered error distribution for the retract
    plt.subplot(gs[2,1])
    plot_filtered_retract_stdev(x_plot,stdev_plot,slice_before,slice_after,
                                epsilon_plot,sigma_plot,ylim_diff_filtered)
    # probability distribution for the retract
    plt.subplot(gs[3,1])
    plot_probability(threshold,x_plot,prob_final,slice_before,slice_after)
    PlotUtilities.label_tom(fig,loc=(-0.15,1.0))
    PlotUtilities.savefig(fig,out_fig,
                          subplots_adjust=dict(hspace=0.2,wspace=0.2))

if __name__ == "__main__":
    run()
