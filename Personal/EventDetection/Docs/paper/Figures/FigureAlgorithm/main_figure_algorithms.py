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
from Research.Personal.EventDetection._2SplineEventDetector import Detector,\
    _no_event

from Research.Personal.EventDetection.Util import Analysis,Plotting

import matplotlib.gridspec as gridspec

# style keywords for this file only
common = dict(rasterized=True)
style_approach = dict(color='k', **common)
style_raw = dict(alpha=0.3,**common)
style_filtered = dict(linewidth=3,**common)
retract_style = dict(color='g',**common)
style_retract_error_dist = dict(linewidth=2,alpha=0.5,**common)
label_s_t = r"$\Sigma_\mathrm{t}$"
label_r_t = r"r$_t$"
lazy_kwargs = dict(loc='lower right',frameon=False)

error_label = "Residual\n(pN)"
filtered_error_label = label_s_t + " (pN)"
force_label = "Force (pN)"

def tick_function():
    tick_kwargs = dict(num_x_major=4,num_y_major=4)
    PlotUtilities.tick_axis_number(**tick_kwargs)


def plot_fec(time_approach,force_approach,interp_approach,
             ylim_force,xlim_approach,label):
    plt.plot(time_approach,force_approach,label=label,
             alpha=0.3,**style_approach)
    plt.plot(time_approach,interp_approach,label="Spline",
             **style_approach)
    PlotUtilities.lazyLabel("",force_label,
                            "Estimating $\epsilon$ and $\sigma$",
                            frameon=False)
    plt.xlim(xlim_approach)
    plt.ylim(ylim_force)

def plot_retract_fec(x_plot,force_plot,slice_before,slice_after,
                     force_filtered_plot,ylim_force):
    Plotting.before_and_after(x_plot,force_plot,slice_before,slice_after,
                              style_raw,label="Retract",label_for_before=False)
    Plotting.before_and_after(x_plot,force_filtered_plot,
                              slice_before,slice_after,style_filtered,
                              label_for_before=False,label="Spline")
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
    plt.plot(time_approach_plot,approach_stdev_plot,alpha=0.3,**style_approach)
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



def plot_epsilon(epsilon_plot,sigma_plot,label="$\epsilon \pm \sigma$"):
    epsilon_style = dict(color='k')
    min_v = epsilon_plot-sigma_plot
    max_v = epsilon_plot+sigma_plot
    plt.axhspan(min_v,max_v,color='c',alpha=0.3,label=label)

def retract_figure(x_plot,force_plot,slice_before,slice_after,
                   force_filtered_plot,ylim_force,force_label,
                   diff_pN,stdev_plot,ylim_diff,ylim_diff_filtered,
                   prob,threshold,epsilon_plot,sigma_plot):
    plt.subplot(4,1,1)
    plot_retract_fec(x_plot,force_plot,slice_before,slice_after,
                     force_filtered_plot,ylim_force)
    # remove the title 
    plt.title("")
    PlotUtilities.xlabel("Time (s)")
    PlotUtilities.ylabel(force_label)
    PlotUtilities.legend(frameon=True,loc='lower right')
    PlotUtilities.x_label_on_top()
    plt.subplot(4,1,2)
    plot_retract_error(x_plot,diff_pN,slice_before,slice_after,stdev_plot,
                       ylim_diff)
    PlotUtilities.ylabel(error_label)
    PlotUtilities.legend(frameon=True,loc='lower right')
    plt.subplot(4,1,3)
    plot_filtered_retract_stdev(x_plot,stdev_plot,slice_before,slice_after,
                                epsilon_plot,sigma_plot,ylim_diff_filtered)
    PlotUtilities.ylabel(filtered_error_label)
    PlotUtilities.legend(frameon=True,loc='upper right')
    plt.subplot(4,1,4)
    plot_probability(threshold,x_plot,prob,slice_before,slice_after)
    PlotUtilities.no_x_label()
    PlotUtilities.xlabel("")

def spline_figures(fec_split,out_base,threshold=1e-2):
    retract = fec_split.retract
    x = retract.Time
    f = retract.Force
    interp = fec_split.retract_spline_interpolator()
    interp_f = interp(x)
    n_points = fec_split.tau_num_points
    delta = _no_event._delta(x,interp_f,n_points)
    approach_kwargs = Detector.make_event_parameters_from_split_fec(fec_split)
    param_obj = _no_event.no_event_parameters(threshold=threshold,
                                              **approach_kwargs)
    # get the delta probability for negative and positive stuff. 
    param_obj.negative_only = False
    delta_prob_all = _no_event._delta_probability(delta,param_obj)
    param_obj.negative_only = True
    delta_prob_neg = _no_event._delta_probability(delta,param_obj)
    y_plot_f = lambda x : x * 1e12
    delta_plot = y_plot_f(delta)
    delta_label = "dF (pN)"
    delta_params = [y_plot_f(param_obj.epsilon),y_plot_f(param_obj.sigma)]
    # get the negtive df parametrs; this is what we actually use
    # (ie: df should be *below* the epsilon and sigm parameters)
    delta_params_n = list(-1 * np.array(delta_params))
    delta_p_label = r"$(\epsilon \pm \sigma)$"
    delta_n_label = r"$-(\epsilon \pm \sigma)$"
    # get the derivative stuff 
    deriv_params = [y_plot_f(param_obj.derivative_epsilon),
                    y_plot_f(param_obj.derivative_sigma)]
    derivative = _no_event._spline_derivative(x,interp)
    deriv_plot = y_plot_f(derivative)
    param_obj.negative_only = True
    deriv_prob = _no_event._derivative_probability(derivative,param_obj)
    deriv_label = r"$\frac{dF}{dt}$ (pN/s)"
    deriv_p_label = r"$\epsilon_{\frac{dF}{dt}} \pm \sigma_{\frac{dF}{dt}}$"
    # get the integral stuff 
    integral = _no_event.local_noise_integral(f,interp_f,n_points)
    integ_plot = y_plot_f(integral)
    integ_prob = _no_event._integral_probability(f,interp_f,n_points,param_obj)
    integ_params = [y_plot_f(param_obj.integral_epsilon),
                    y_plot_f(param_obj.integral_sigma)]
    integ_label = r"$\int r_t dt$ (pN $\cdot$ s)"
    integ_p_label = r"$\epsilon_{\int} \pm \sigma_{\int}$"
    # plot all the values and their distributions
    plot_tuples = [ \
        ([delta_plot,delta_prob_all,delta_label,delta_p_label] +delta_params),
        ([delta_plot,delta_prob_neg,delta_label,delta_n_label] +delta_params_n),
        ([deriv_plot,deriv_prob    ,deriv_label,deriv_p_label] +deriv_params),
        ([integ_plot,integ_prob    ,integ_label,integ_p_label] +integ_params),
                    ]
    style_plotting = dict(color='k')
    for i,(val,prob,y_label,eps_label,eps,sigma) in enumerate(plot_tuples):
        fig = PlotUtilities.figure((8,6))
        plt.subplot(2,1,1)
        plt.plot(x,val,**style_plotting)
        plot_epsilon(eps,sigma,label=eps_label)
        PlotUtilities.lazyLabel("Time (s)",y_label,"")
        PlotUtilities.x_label_on_top()
        plt.subplot(2,1,2)
        plt.semilogy(x,prob,**style_plotting)
        PlotUtilities.lazyLabel("","Probability","")
        PlotUtilities.no_x_label()
        PlotUtilities.savefig(fig,out_base + "spline{:d}.pdf".format(i))

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
            stdevs_epsilon_and_sigma(slice_fit_approach=slice_fit_approach)
    # get the retract 
    interpolator = fec_split.retract_spline_interpolator()
    force_filtered = interpolator(time)
    n_points = fec_split.tau_num_points
    diff = force-force_filtered
    diff_pN = y_plot_f(diff)
    stdev,_,_ = Analysis.stdevs_epsilon_sigma(force,force_filtered,
                                              n_points)
    _,epsilon,sigma = fec_split.\
        stdevs_epsilon_and_sigma(n_points=n_points,
                                slice_fit_approach=slice_fit_approach)
    # get the probability
    threshold = 1e-2
    approach_kwargs = Detector.make_event_parameters_from_split_fec(fec_split)
    # get the prediction info
    _,predict_info = Detector._predict_full(example,threshold=threshold)
    # get the final masks
    mask_final = predict_info.condition_results[-1]
    bool_final = np.zeros(time.size)
    bool_final[mask_final] = 1
    prob_final = predict_info.probabilities[-1]
    prob_limits =  [min(prob_final/2),2]
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
    n_rows = 3
    n_cols = 2
    subplots_adjust_pres = dict(hspace=0.1)
    """
    Make the spline probability figures
    """
    spline_figures(fec_split,out_base="./explaining_")
    """
    make just the retract figures
    """
    retract_figsize = (6,8)
    # make retract figures showing how the probability is built up
    extra_kwargs = [dict(valid_delta=False,valid_integral=False,
                         valid_derivative=False),
                    dict(valid_delta=False,valid_integral=False),
                    dict(valid_delta=False),
                    dict()]
    for i,extra in enumerate(extra_kwargs):
        full_dict = dict(threshold=threshold,**approach_kwargs)
        full_dict = dict(full_dict,negative_only=True,mask_is_conditional=False,
                         **extra)
        obj_tmp = _no_event.no_event_parameters(**full_dict)
        prob_tmp,_= _no_event.\
        _no_event_probability(time,interp=interpolator,
                              y=force,n_points=n_points,
                              no_event_parameters_object=obj_tmp)
        fig = PlotUtilities.figure(retract_figsize)
        retract_figure(x_plot,force_plot,slice_before,slice_after,
                       force_filtered_plot,ylim_force,force_label,
                       diff_pN,stdev_plot,ylim_diff,ylim_diff_filtered,
                       prob_tmp,threshold,epsilon_plot,sigma_plot)
        plt.ylim(prob_limits)
        out_replace = "_retract{:d}.pdf".format(i)
        PlotUtilities.savefig(fig,out_fig.replace(".pdf",out_replace),
                              subplots_adjust=subplots_adjust_pres)
    f_adhesions = Detector.adhesion_mask_function_for_split_fec
    f_delta = Detector.delta_mask_function
    # make the 'final result' figure, showing how the masks are used
    kwargs_masks = [dict(f_refs=[f_adhesions],mask_is_conditional=False),
                    dict(f_refs=[f_adhesions,f_delta],
                                 mask_is_conditional=False)]
    offset = len(extra_kwargs)
    for i,kwargs_tmp in enumerate(kwargs_masks): 
         _,predict_info_tmp = Detector._predict_full(example,
                                                     threshold=threshold,
                                                     **kwargs_tmp)
         prob_tmp = predict_info_tmp.probabilities[-1]
         fig = PlotUtilities.figure(retract_figsize)
         retract_figure(x_plot,force_plot,slice_before,slice_after,
                        force_filtered_plot,ylim_force,force_label,
                        diff_pN,stdev_plot,ylim_diff,ylim_diff_filtered,
                        prob_tmp,threshold,epsilon_plot,sigma_plot)
         plt.ylim(prob_limits)         
         out_replace = "_retract{:d}.pdf".format(i+offset)
         PlotUtilities.savefig(fig,out_fig.replace(".pdf",out_replace),
                               subplots_adjust=subplots_adjust_pres)
    """
    make just the approach figure
    """
    fig = PlotUtilities.figure((6,8))
    plt.subplot(3,1,1)
    plot_fec(time_approach,force_approach,interp_approach,
             ylim_force,xlim_approach,"Approach")
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
    PlotUtilities.savefig(fig,out_fig.replace(".pdf","_approach.pdf"),
                          subplots_adjust=subplots_adjust_pres)
    """
    make the plot for the paper
    """
    fig = PlotUtilities.figure((7,5))
    gs = gridspec.GridSpec(4, 2)
    plt.subplot(gs[0,0])
    plot_fec(time_approach,force_approach,interp_approach,
             ylim_force,xlim_approach,"Approach")
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
    PlotUtilities.label_tom(fig,loc=(-0.1,1.05))
    PlotUtilities.save_twice(fig,out_fig + ".png",out_fig+".svg",
                             subplots_adjust=dict(hspace=0.3,wspace=0.15))


if __name__ == "__main__":
    run()
