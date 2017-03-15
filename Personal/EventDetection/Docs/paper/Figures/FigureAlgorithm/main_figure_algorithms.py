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

def fmt(ax):
    ax.set_xlim([-0.1,2.5])
    ax.set_ylim([-50,40])

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
    # get the approach 
    approach = fec_split.approach
    approach_time = approach.Time
    time_approach = x_plot_f(approach_time)
    force_approach = y_plot_f(approach.Force)
    interpolator_approach = fec_split.approach_spline_interpolator()
    interp_approach = y_plot_f(interpolator_approach(approach_time))
    # get the retract 
    interpolator = fec_split.retract_spline_interpolator()
    force_filtered = interpolator(time)
    n_points = fec_split.tau_num_points
    diff = force-force_filtered
    diff_pN = y_plot_f(diff)
    stdev,_,_ = Analysis.stdevs_epsilon_sigma(force,force_filtered,
                                              n_points)
    # get the epsilon and sigmas
    n = fec_split.approach.Force.size
    offset_approach =  n-fec_split.get_predicted_approach_surface_index()
    slice_fit_approach = slice(offset_approach,-offset_approach,1)
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
    force_plot = y_plot_f(force)
    force_filtered_plot = y_plot_f(force_filtered)
    epsilon_plot,sigma_plot = y_plot_f(epsilon),y_plot_f(sigma)
    stdev_plot = y_plot_f(stdev)
    before,after = Analysis.get_before_and_after_and_zoom_of_slice(fec_split)
    slice_before,slice_after = before[0],after[0]
    fig = PlotUtilities.figure((16,8))
    n_rows = 3
    n_cols = 2
    epsilon_style = dict(color='b')
    plt.subplot(n_rows,n_cols,1)
    style_approach = dict(color='k')
    plt.plot(time_approach,force_approach,label="Raw Force (approach)",
             alpha=0.3,**style_approach)
    plt.plot(time_approach,interp_approach,label="Spline ($g^{*}_t$)))",
             **style_approach)
    PlotUtilities.lazyLabel("","Force [pN]","",frameon=True)
    plt.subplot(n_rows,n_cols,2)
    style_raw = dict(alpha=0.3)
    style_filtered = dict(linewidth=3)
    Plotting.before_and_after(x_plot,force_plot,slice_before,slice_after,
                              style_raw,label="Raw Force (retract)")
    Plotting.before_and_after(x_plot,force_filtered_plot,
                              slice_before,slice_after,style_filtered,
                              label="Spline ($g^{*}_t$)))")
    PlotUtilities.lazyLabel("","","",**lazy_kwargs)
    plt.subplot(n_rows,n_cols,3)
    plt.plot(x_plot,diff_pN,alpha=0.3,**epsilon_style)
    plt.plot(x_plot,stdev_plot,**epsilon_style)
    PlotUtilities.lazyLabel("","R$_\mathrm{t}$ [pN]","",frameon=True)
    plt.subplot(n_rows,n_cols,4)
    plt.plot(x_plot,stdev_plot)
    plt.axhline(epsilon_plot,label="$\epsilon$")
    plt.axhline(epsilon_plot+sigma_plot,linestyle='--',
                label="$\epsilon \pm \sigma$",**epsilon_style)
    plt.axhline(epsilon_plot-sigma_plot,linestyle='--')
    PlotUtilities.lazyLabel("","R$_\mathrm{t}$^{*}[pN]","",frameon=True)
    plt.subplot(n_rows,n_cols,6)
    plt.plot(x_plot,prob,alpha=0.3,color='k',label="No-event")
    plt.axhline(threshold,linewidth=2,color='b',linestyle='--',
                label="threshold")
    plt.plot(x_plot,prob_final,label="Masked no-event")
    plt.yscale('log')
    PlotUtilities.lazyLabel("Time (s)","Probability","",**lazy_kwargs)
    PlotUtilities.label_tom(fig,loc=(-1.13,1.05))
    PlotUtilities.savefig(fig,out_fig)
    

if __name__ == "__main__":
    run()
