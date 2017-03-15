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
from Research.Personal.EventDetection._2SplineEventDetector import Detector

from Research.Personal.EventDetection.Util import Analysis,Plotting

import matplotlib.gridspec as gridspec

def fmt(ax):
    ax.set_xlim([-0.1,2.5])
    ax.set_ylim([-50,40])

def plot_epsilon(epsilon_plot,sigma_plot):
    epsilon_style = dict(color='b')
    plt.axhline(epsilon_plot,label="$\epsilon$")
    plt.axhline(epsilon_plot+sigma_plot,linestyle='--',
                label="$\epsilon \pm \sigma$",**epsilon_style)
    plt.axhline(epsilon_plot-sigma_plot,linestyle='--')

def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    out_fig = "algorithm.pdf"
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
    lazy_kwargs = dict(loc='lower right',frameon=True)
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
    fig = PlotUtilities.figure((16,12))
    n_rows = 3
    n_cols = 2
    tick_kwargs = dict(num_x_major=4,num_y_major=4)
    tick_function = lambda: PlotUtilities.tick_axis_number(**tick_kwargs)
    gs = gridspec.GridSpec(4, 2)
    plt.subplot(gs[0,0])
    style_approach = dict(color='k')
    plt.plot(time_approach,force_approach,label="Raw Force (approach)",
             alpha=0.3,**style_approach)
    plt.plot(time_approach,interp_approach,label="Spline ($g^{*}_t$)))",
             **style_approach)
    PlotUtilities.lazyLabel("","Force [pN]","",frameon=True)
    plt.xlim(xlim_approach)
    tick_function()
    PlotUtilities.no_x_ticks()
    plt.subplot(gs[0,1])
    style_raw = dict(alpha=0.3)
    style_filtered = dict(linewidth=3)
    Plotting.before_and_after(x_plot,force_plot,slice_before,slice_after,
                              style_raw,label="Raw Force (retract)")
    Plotting.before_and_after(x_plot,force_filtered_plot,
                              slice_before,slice_after,style_filtered,
                              label="Spline ($g^{*}_t$)))")
    PlotUtilities.lazyLabel("","","",**lazy_kwargs)
    tick_function()
    PlotUtilities.no_x_ticks()
    retract_style = dict(color='g')
    # plot the 'raw' error distribution for the approach
    plt.subplot(gs[1,0])
    plt.plot(time_approach_plot,approach_diff[slice_fit_approach],alpha=0.3,
             **style_approach)
    plt.plot(time_approach_plot,approach_stdev_plot,**style_approach)
    PlotUtilities.lazyLabel("","R$_\mathrm{t}$ [pN]","",frameon=True)
    plt.xlim(xlim_approach)
    PlotUtilities.no_x_ticks()
    tick_function()
    plt.ylim(*ylim_diff)
    # plot the 'raw error distribution for the retract
    plt.subplot(gs[1,1])
    plt.plot(x_plot,diff_pN,alpha=0.3,**retract_style)
    plt.plot(x_plot,stdev_plot,**retract_style)
    PlotUtilities.lazyLabel("","","")
    PlotUtilities.no_x_ticks()
    plt.ylim(*ylim_diff)
    tick_function()
    # filtered error distribution for the approach
    plt.subplot(gs[2,0])
    plt.plot(time_approach_plot,approach_stdev_plot,**style_approach)
    plot_epsilon(epsilon_plot,sigma_plot)
    plt.xlim(xlim_approach)
    plt.ylim(*ylim_diff_filtered)
    tick_function()
    PlotUtilities.lazyLabel("Time (s)","R$_{\mathrm{t}}$ [pN]","",
                            frameon=True,loc='upper right')
    # filtered error distribution for the retract
    plt.subplot(gs[2,1])
    plt.plot(x_plot,stdev_plot,**retract_style)
    plot_epsilon(epsilon_plot,sigma_plot)
    PlotUtilities.lazyLabel("","","",frameon=True,
                            loc="upper right")
    PlotUtilities.no_x_ticks()
    plt.ylim(*ylim_diff_filtered)
    tick_function()
    # probability distribution for the retract
    plt.subplot(gs[3,1])
    plt.yscale('log')
    plt.plot(x_plot,prob,alpha=0.3,label="No-event",**retract_style)
    plt.axhline(threshold,linewidth=3,linestyle='--',
                label="threshold",color='k')
    plt.plot(x_plot,prob_final,label="Masked no-event",**retract_style)
    PlotUtilities.lazyLabel("Time (s)","Probability","",**lazy_kwargs)
    tick_function()
    PlotUtilities.label_tom(fig,loc=(-1.10,1.0))
    PlotUtilities.savefig(fig,out_fig,
                          subplots_adjust=dict(hspace=0.2,wspace=0.2))
    

if __name__ == "__main__":
    run()
