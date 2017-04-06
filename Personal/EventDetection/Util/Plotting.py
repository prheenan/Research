# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import interpolate
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from GeneralUtil.python import PlotUtilities
from scipy.stats import norm
from Research.Personal.EventDetection.Util import Analysis,Learning,Offline
import matplotlib.gridspec as gridspec

style_train = dict(color='r',marker='o',linestyle='--',label="Training") 
style_valid = dict(color='g',marker='v',linestyle='-',label="Validation")
color_pred_def = 'b'
color_true_def = 'g'

def _style_true(color_true=color_true_def):
    return dict(color=color_true,alpha=0.9)

def _style_pred(color_pred=color_pred_def):
    return dict(color=color_pred,alpha=0.4)

_fec_event_colors = ['k','r','b']

def algorithm_colors():
    return ['b','k','r']

def algorithm_markers():
    return ['s','o','v']

def algorithm_linestyles():
    return ['--','-.','-']    

def algorithm_title_dict():
    """
    returns a dicitonary mapping <description of algorithm> -> plot title
    """
    dict_v = {"no event":"FEATHER",
              "open fovea":"OpenFovea",
              "wavelet transform":"Scientific Python"}
    return dict_v

def true_hatch():
    return "//"

def plot_autocorrelation_log(x,*args):
    """
    plots the autocorrelation function and fit
    
    Args:
        x: corrlation abscissa
        *args: output of auto_correlation_tau
    Returns:
        nothing, plots the autocorrelation log 
    """
    tau,coeffs,auto = args
    auto_norm_predicted = np.exp(np.polyval(coeffs,x=x))
    auto_norm,statistical_norm,log_norm,fit_idx_max= \
        Analysis.auto_correlation_helper(auto)
    tol = 1e-6
    x_plot = x- min(x)
    # auto_norm goes from 0 to 1
    style_fit = dict(color='r',linestyle='--',linewidth=2)
    style_auto = dict(color='k',alpha=0.3,marker='.')
    xlim_zoomed = [0,3*tau]
    highlight_fit_range = lambda y: \
        plt.plot(x_plot[:fit_idx_max],y[:fit_idx_max],color='r')
    plt.subplot(3,1,1)
    plt.plot(x_plot,auto_norm,**style_auto)
    plt.plot(x_plot,auto_norm_predicted,**style_fit)
    highlight_fit_range(auto_norm)
    PlotUtilities.lazyLabel("","autocorrelation (au)","")
    plt.subplot(3,1,2)
    plt.plot(x_plot,auto_norm,**style_auto)
    plt.plot(x_plot,auto_norm_predicted,**style_fit)
    highlight_fit_range(auto_norm)
    plt.xlim(xlim_zoomed)
    PlotUtilities.lazyLabel("","","")
    plt.subplot(3,1,3)
    plt.plot(x_plot,log_norm,**style_auto)
    plt.plot(x_plot,np.polyval(coeffs,x=x),**style_fit)
    plt.xlim(xlim_zoomed)
    plt.ylim([np.percentilepercentile(log_norm,0.5),max(log_norm)])
    PlotUtilities.lazyLabel("time","log of normalized autocorrelation","")

         
def plot_autocorrelation(example_split): 
    """
    Given an already-split force extension curve, plots the autocorrelation
    information related to it.

    Args:
        example_splut: split_force_extesion object, already zeroed
    Returns:
        nothing, sets up a plott..
    """
    retract = example_split.retract
    x,separation,f = retract.Time,retract.Separation,retract.Force
    tau = example_split.tau
    # XXX only look at after the nominal zero point?
    # get an interpolator for the retract force and separation
    force_interpolator = Analysis.spline_interpolator(tau,x,f)
    separation_interpolate = Analysis.spline_interpolator(tau,x,separation)
    tau,auto_coeffs,auto_correlation = Analysis.auto_correlation_tau(x,f)
    # plot everything
    plot_autocorrelation_log(x, tau,auto_coeffs,auto_correlation)

def highlight_events(idx_events,x,y,label=None,**kwargs):
    """
    highlights x and y at the indices given by idx_events:
    
    Args:
        idx_events: where the events are (a list of slices or integers
        x,y: what to plot
        label: how to label it
        **kwargs: passed to highlight_events
    """
    for i in idx_events:
        # only label the first
        label_tmp = label if i == 0 else None
        plt.plot(x[i],y[i],label=label_tmp,**kwargs)

        
def plot_prediction_info(ex,info,xlabel="Time",
                         x_func=lambda ex: ex.Time -ex.Time[0]):   
    """
    plots prediction information for the no-hypothesis predictor

    Args:
        ex: instance of (zeroed, etc) TimeSepForce instance
    
    Returns:
        This is a description of what is returned.
    """
    thresh = info.threshold
    retract = ex.retract
    event_slices = ex.get_retract_event_idx()
    time,separation,force = retract.Time,retract.Separation,retract.Force
    event_slices_predicted = info.event_slices
    event_idx_start = [e.start for e in event_slices_predicted]
    event_idx_end = [e.stop for e in event_slices_predicted]
    event_idx = info.event_idx
    mask = info.mask
    # get the interpolated derivative
    slice_v = info.slice_fit
    time_slice = time[slice_v]
    interpolator = ex.retract_spline_interpolator(slice_v)
    interp_first_deriv = interpolator.derivative(1)(time_slice)
    interpolated_force = interpolator(time_slice)
    tau = ex.tau
    stdevs = info.local_stdev
    # plot everything
    style_events = dict(color='r',label="True events")
    x = x_func(retract)
    surface_index = ex.get_predicted_retract_surface_index()
    min_x,max_x = min(x),max(x)
    ylabel = "Force (pN)"
    x_range = max_x - min_x
    fudge = x_range * 0.05
    x_limits = [min_x - fudge,max_x + fudge]
    force_plot = force * 1e12
    interpolated_force_plot = interpolated_force*1e12
    time_interpolated_plot = time_slice - min(retract.Time)
    # get the informaiton relevant to the CDF
    original_cdf = info.probabilities[0]
    cdf = original_cdf
    boolean_mask = np.zeros_like(cdf)
    boolean_mask[mask] = 1
    masked_cdf = cdf.copy()
    masked_cdf *= boolean_mask
    n_rows = 3
    n_cols = 1
    lazy_kwargs = dict(frameon=False,loc="best")
    plt.subplot(n_rows,n_cols,1)
    plt.plot(x,force_plot,color='k',alpha=0.3,label="data")
    plt.plot(time_interpolated_plot,interpolated_force_plot,color='b',
             linewidth=2,label="2-spline")
    plt.axvline(x[surface_index],label="Surface\n(pred)")
    highlight_events(event_slices,x,force_plot,**style_events)
    PlotUtilities.lazyLabel("",ylabel,"",**lazy_kwargs)
    plt.xlim(x_limits)
    plt.subplot(n_rows,n_cols,2)
    # plot the autocorrelation time along the plot
    min_x_auto = min(x) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+tau]
    styles = [dict(color='k',linestyle='-',alpha=0.3),
              dict(color='b',linestyle='-.',alpha=0.3),
              dict(color='y',linestyle=':',alpha=1),
              dict(color='g',linestyle='-',alpha=0.7),
              dict(color='m',linestyle='--'),]
    for i,c in enumerate(info.probabilities):
        sty = styles[i % len(styles)]
        plt.semilogy(x,c,label="cdf{:d}".format(i),**sty)
    min_cdf = min([min(c) for c in info.probabilities])
    plt.axhline(thresh,color='k',linestyle='--',label="tsh")
    mask_boolean = np.zeros(x.size)
    mask_boolean[mask] = 1
    PlotUtilities.lazyLabel("","No-Event CDF ","",loc='upper right',
                            frameon=False)
    plt.xlim(x_limits)
    plt.ylim([min_cdf/2,3])
    mask_styles = styles
    tol = 1.5 + 1e-6
    for i,c in enumerate(info.condition_results):
        bool = np.zeros_like(x)
        bool[c] = 1
        style = mask_styles[i % len(mask_styles)]
        plt.semilogy(x,bool+tol,label="mask {:d}".format(i),**style)
    plt.plot(x,boolean_mask,color='k',linestyle='-',label="final mask")
    plt.xlim(x_limits)    
    plt.subplot(n_rows,n_cols,3)
    # XXX check mask has at least one...
    plt.plot(x,force_plot,color='k',alpha=0.3,label="data")
    plt.plot(x[mask],force_plot[mask],color='k',alpha=0.8,
             label="Event region")
    for i,(fwd,rev,event) in \
        enumerate(zip(event_idx_end,event_idx_start,event_idx)):
        try:
            plt.axvline(x[fwd],linestyle='--',color='r')
            plt.axvline(x[rev],linestyle='-.',color='g')
            label = "predicted" if (i ==0) else ""            
            plt.plot(x[event],force_plot[event],marker='o',color='b',alpha=0.3,
                    label=label)
        except IndexError as e:
            print(e)
    plt.xlim(x_limits)
    PlotUtilities.lazyLabel(xlabel,ylabel,"",**lazy_kwargs)

def plot_classification(split_object,scoring_object):
    """
    plots the classification of a single object

    Args:
         split_object: whatever was classified
         scoring_object: Scoring.score object, used on split_object
    Returns:
         Nothing
    """
    retract = split_object.retract
    time,force = retract.Time,retract.Force
    true_idx,pred_idx =  scoring_object.idx_true,scoring_object.idx_predicted
    boolean_true,boolean_predicted = scoring_object.get_boolean_arrays(time)
    plt.subplot(2,1,1)
    plt.plot(time,force,color='k',alpha=0.3)
    highlight_events(true_idx,time,force,color='r',linestyle='-',
                     label="true events")
    highlight_events(pred_idx,time,force,color='b',linestyle='None',
                     marker='.',label="predicted events")
    PlotUtilities.lazyLabel("","Force(pN)","")
    plt.subplot(2,1,2)
    plt.plot(time,boolean_true,linewidth=4,label="True events",color='b')
    plt.plot(time,boolean_predicted,label="Predicted",color='r',alpha=0.3)
    PlotUtilities.lazyLabel("Time (s)","Force(pN)","")

def debug_plot_event(x,y,fit_x,fit_y,x_event,y_event,pred,idx_above_predicted):
    """
    Used  insides of Detector.event_by_loadin_rate to plot what is happening
    
    Args:
        all: internal Detector.event_by_loading_rate, see that function
    Returns:
        nothing
    """          
    x_plot = lambda x_tmp: x_tmp - min(x)
    y_plot = lambda y_tmp : y_tmp * 1e12
    plt.subplot(2,1,1)
    plt.plot(x_plot(x),y_plot(y),alpha=0.3,color='k',label="raw")
    plt.plot(x_plot(fit_x),y_plot(fit_y),label="fit")
    PlotUtilities.lazyLabel("","Force (pN)","",frameon=False)
    plt.subplot(2,1,2)
    plt.plot(x_plot(fit_x),y_plot(fit_y))
    plt.plot(x_plot(x_event),y_plot(y_event),color='r',alpha=0.3,label="event")
    plt.plot(x_plot(x_event),y_plot(pred),label="prediction")
    if (len(idx_above_predicted) > 0):
        plt.axvline(x_plot(x[idx_above_predicted[-1]]))
    PlotUtilities.lazyLabel("Time","Force (pN)","",frameon=False,
                            loc="upper left")

def debug_plot_adhesion_info(split_fec,min_idx,boolean_ret):
    """
    Used inside of Detector.adhesion_mask to tell wtf is happening
    
    Args:
        all: internal Detector.adhesion_mask, see that function
    Returns:
        nothing
    """      
    surface_index = split_fec.get_predicted_retract_surface_index()    
    x = split_fec.retract.Time
    force = split_fec.retract.Force
    plt.subplot(2,1,1)
    plt.plot(x,force*1e12,color='k',alpha=0.3)
    plt.axvline(x[min_idx],label="minimum")
    plt.axvline(x[surface_index],label="surface",linestyle='--')
    PlotUtilities.lazyLabel("","Force","",loc="upper right",
                            frameon=True)     
    plt.subplot(2,1,2)
    plt.plot(x,boolean_ret,color='k',linestyle='--')
    PlotUtilities.lazyLabel("Time","mask","",loc="upper right",
                            frameon=True)   

def _plot_rupture_objects(to_plot,**kwargs):
    """
    given rupture objects, plots rupture force in pN (given N) vs
    log of loading rate in pN/s (given N/s)

    Args:
         to_plot: list of rupture obects to plot with the same style
         **kwargs: passed to plt.semilogx
    Returns:
         Nothing
    """
    rupture_forces_pN,loading_rate_pN_per_s = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(to_plot)
    plt.semilogx(loading_rate_pN_per_s,rupture_forces_pN,**kwargs)

def plot_true_and_predicted_ruptures(true,predicted,title="",
                                     loc='upper left',
                                     style_predicted=None,style_true=None):
    """
    given rupture objects, plots the true and predicted values of rupture
    force verus loading rate

    Args:
        true / predicted: list of true and predicted ruptures. dont have to 
        match up, but should be all from the same FECs

        title: of the plot
        label_true: for the legend marker on the true objects
        style_predicted: what to make the predicted ones look like. if None,
        defsults to little blue x's.
    Returns:
         Nothing
    """
    if (style_predicted is None):
        style_predicted = dict(label="predicted",linewidth=2,
                               **_style_pred_def('k'))
    if (style_true is None):
        style_true = dict(label="true",**_style_true_def('g'))
    marker_size = 4
    style_true_marker = dict(**style_true)
    style_true_marker['alpha'] = 0.5
    _plot_rupture_objects(predicted,marker='x',linewidth=3,linestyle='None',
                          markersize=marker_size,color=style_predicted['color'],
                          label=style_predicted['label'])
    _plot_rupture_objects(true,marker='o',linewidth=4,markersize=marker_size,
                          markerfacecolor="None",markeredgecolor='g',
                          linestyle='None',**style_true_marker)
    PlotUtilities.lazyLabel("Loading Rate (pN/s)","Rupture Force (pN)",
                            title,frameon=True,
                            legend_kwargs=dict(numpoints=1,markerscale=2),
                            loc=loc)
    PlotUtilities.set_legend_kwargs()

def debug_plot_signal_mask(x,force,gt_condition,x_sliced,interp_f,
                           boolean_array,no_event_cond,value_cond,
                           boolean_ret,probability_updated,probability,
                           threshold):
    xlim = plt.xlim(min(x),max(x))
    plt.subplot(4,1,1)
    valid_idx = np.where(np.logical_not(gt_condition))
    invalid_idx = np.where(gt_condition)
    plt.plot(x[invalid_idx],force[invalid_idx],color='k',alpha=0.3)
    plt.plot(x[valid_idx],force[valid_idx],color='g')    
    plt.plot(x_sliced,interp_f,color='b')
    plt.xlim(xlim)
    plt.subplot(4,1,2)
    plt.plot(x,boolean_array+2.1,label="orig")
    plt.plot(x_sliced,no_event_cond+1.1)
    plt.plot(x_sliced,value_cond)
    plt.plot(x,gt_condition-1.1,label="concat")
    plt.plot(x,boolean_ret-2.1,linestyle='--',label="f")
    plt.legend(loc='upper left')
    plt.xlim(xlim)
    plt.subplot(4,1,3)
    plt.plot(x,boolean_array+1.1)
    plt.plot(x,boolean_ret+2.1,linestyle='--')
    plt.xlim(xlim)
    plt.subplot(4,1,4)
    plt.semilogy(x,probability_updated,linestyle='--')
    plt.semilogy(x,probability)
    plt.xlim(xlim)
    plt.axhline(threshold)
    plt.xlim(xlim)
    
    
def debugging_plots(id_string,example_split,info,plot_auto=False):
    """
    Plots the autocorrelation and prediction information

    Args:
        id_string: 
        example_split: force extension curve to debug
        info: return from predict_helper
    Returns:
        nothing, splits plots out like id_string
    """
    out_file_path =  id_string
    if (plot_auto):
        fig = PlotUtilities.figure(figsize=(8,20))
        plot_autocorrelation(example_split)
        PlotUtilities.savefig(fig,out_file_path + "auto.png")   
    # XXX fix threshhold
    fig = PlotUtilities.figure(figsize=(8,12))    
    plot_prediction_info(example_split,info)
    PlotUtilities.savefig(fig,out_file_path + "info.png")

def plot_ruptures_of_scores(scores):
    """
    Given a list of score objects, plots the rupture distributions

    Args:
         scorer objects
    Returns:
         nothing
    """
    rupture_true = [r for s in scores for r in s.ruptures_true ]
    rupture_predicted = [r for s in scores for r in s.ruptures_predicted ]
    plot_true_and_predicted_ruptures(rupture_true,rupture_predicted)

def cross_validation_distance_metric(x_values,train_scores,valid_scores,
                                     to_true):
    """
    Plots the cross validation training and validation distance metric

    Args:
        x_values: x values
        <train/valid>_scores: the training and validation scores used
        to_true: distance metric plotted is *from* predicted *to* true 
        (ie: something like recall) if true, otherwise vice versa

    Returns:
        nothing
    """
    # get the metrics and errors by parameters
    x_train,train_dist,train_dist_std = \
        Learning.median_dist_metric(x_values,train_scores,to_true=to_true)
    x_valid,valid_dist,valid_dist_std = \
        Learning.median_dist_metric(x_values,valid_scores,to_true=to_true)
    y_plot = lambda y: y * 1e9
    train_dist_plot,train_error_plot = y_plot(train_dist),y_plot(train_dist_std)
    valid_dist_plot,valid_error_plot = y_plot(valid_dist),y_plot(valid_dist_std)
    plt.errorbar(x=x_train,y=train_dist_plot,yerr=train_error_plot,
                 **style_train)
    plt.errorbar(x=x_valid,y=valid_dist_plot,yerr=valid_error_plot,
                 **style_valid)
    plt.xscale('log')
    plt.yscale('log')
    PlotUtilities.lazyLabel("Tuning Parameter","Median event distance (nm)","",
                            frameon=False)
    
def get_train_test_n_off_and_error(x_values,train_scores,valid_scores):
    x_train,train_dist,train_dist_std = \
        Learning.number_events_off_per_param(x_values,train_scores)
    x_valid,valid_dist,valid_dist_std = \
        Learning.number_events_off_per_param(x_values,valid_scores)
    return x_train,train_dist,train_dist_std,x_valid,valid_dist,valid_dist_std

def plot_num_events_off(x_values,train_scores,valid_scores,ylim=None):
    """
    Plots the number of extra or missing events (irrespective of where they are0

    Args:
        cache_directory: where to save the plots
        learner: learning_curve instance to use
    Returns:
        nothing
    """
    x_train,train_dist,train_dist_std,x_valid,valid_dist,valid_dist_std = \
        get_train_test_n_off_and_error(x_values,train_scores,valid_scores)
    _plot_num_events_off(x_train,train_dist,train_dist_std,
                         x_valid,valid_dist,valid_dist_std)
    if ylim is None:
        ylim = [1e-2,max(train_dist+train_dist_std)]
    plt.ylim(ylim)

def _plot_num_events_off(x_train,train_dist,train_error,
                          x_valid,valid_dist,valid_error,ylim=None,
                         xlabel=None,ylabel=None,lazy_kwargs=dict()):
    """
    see plot_num_events off

    Args;
        x_<y> : the x values for data type y (valid or train)
        <y>_dist: the distribution of values of y (ie: graphical y values)
        <y>_error: the error of the graphical y values
        others: plotting
    """
    if (xlabel is None):
        xlabel = "Tuning parameter"
    if (ylabel is None):
        ylabel = "Relative number of missing or incorrect events"
    plt.errorbar(x=x_train,y=train_dist,yerr=train_error,**style_train)
    plt.errorbar(x=x_valid,y=valid_dist,yerr=valid_error,**style_valid)
    PlotUtilities.lazyLabel(xlabel,ylabel,"",**lazy_kwargs)
    plt.xscale('log')    
    plt.yscale('log')    

def distance_distribution_plot(learner,box_kwargs=None,**kwargs):
    """
    plots the distribution of distances to/from predicted events from/to
    actual events, dependning on kwargs
    
    Args:
        learner: the learner object to use
        kwargs: passed to event_distance_distribution (ie: to_true=T/F)
    """
    train_scores = learner._scores_by_params(train=True)
    valid_scores = learner._scores_by_params(train=False)
    if (box_kwargs is None):
        box_kwargs = dict(whis=[5,95])
    name = learner.description.lower()
    x_values = learner.param_values()
    train_dist = Learning.event_distance_distribution(train_scores,**kwargs)
    valid_dist = Learning.event_distance_distribution(valid_scores,**kwargs)
    dist_plot = lambda x: [v * 1e9 for v in x]
    train_plot = dist_plot(train_dist)
    valid_plot = dist_plot(valid_dist)
    plt.boxplot(x=train_plot,**box_kwargs)
    plt.boxplot(x=valid_plot,**box_kwargs)
    plt.gca().set_yscale('log')
    PlotUtilities.lazyLabel("Tuning parameter","Distance Distribution (nm)",
                            "Event distributions for {:s}".format(name),
                            frameon=False)


def _histogram_true_style(color_true=color_true_def,label="True"):
    style_true = dict(color=color_true,label=label,edgecolor=color_true,
                      histtype='stepfilled',fill=True,hatch= true_hatch(),
                      facecolor=color_true,alpha=0.4)
    return style_true

def _histogram_predicted_style(color_pred=color_pred_def,label="Predicted"):
    style_pred = dict(color=color_pred,label=label,fill=False,
                      histtype='step',alpha=1,linewidth=3)
    return style_pred

def event_error_kwargs(metric,color_pred='b',color_true='g',n_bins = 50,
                       xlabel="Relative Error (x$_\mathrm{k}$)",
                       distance_limits=None,clip_limits=False,q=None):
    """
    Args:
        see Plotting.histogram_event_distribution
    Returns:
        a dict with the event error kwargs, see 
        Plotting.histogram_event_distribution
    """
    name = metric.name.lower()
    if (name == "no event"):
        loc = "upper right"
    else:
        loc = "upper left"
    label_pred = r"d$_{\mathrm{p}\rightarrow\mathrm{t}}$"
    label_true = r"d$_{\mathrm{t}\rightarrow\mathrm{p}}$"
    style_pred = _histogram_predicted_style(color_pred=color_pred,
                                            label=label_pred)
    style_true = _histogram_true_style(color_true=color_true,label=label_true)

    to_true,to_pred = metric.to_true_and_pred_distances()
    limit = metric.distance_limit(relative=True)
    log_limit = np.log10(limit)
    max_x_true,max_x_pred =  metric.max_x_distances_true_pred()
    bins = np.logspace(*log_limit,num=n_bins)
    if (distance_limits is None):
        distance_limits = limit
    return dict(to_true=to_true,to_pred=to_pred,distance_limits=distance_limits,
                bins=bins,style_true=style_true,style_pred=style_pred,loc=loc,
                xlabel=xlabel,max_x_true=max_x_true,max_x_pred=max_x_pred,
                q=q)


def histogram_event_distribution(to_true,to_pred,distance_limits,bins,
                                 style_true,style_pred,max_x_true,max_x_pred,
                                 xlabel="Distance [m]",loc='best',q=None,
                                 q_label=None):
    """
    plots the distribution of distances from true/predicted to counterparts

    Args:
        to_<y> : list of distances to y from its counterpart
        distance_limits: x limits
        bins: fed to plt.hist
        style_<x>:  the histogram style for x
        xlabel: label for the x axis 
    """
    # plot the distance scores; color by i in 'i->j' (ie: true->predicted
    # is colored by true
    rel_pred = to_pred/max_x_true
    rel_true = to_true/max_x_pred
    if (to_pred.size > 0):
        plt.hist(rel_pred,
                 log=True,bins=bins,**style_true)
    if (to_true.size > 0):
        plt.hist(rel_true,log=True,bins=bins,**style_pred)
    if (q is not None):
        cat = np.concatenate([rel_pred,rel_true])
        q_num = np.percentile(cat,q)
        if (q_label is None):
            q_label = ((r"P$_{" + "{:d}".format(q) + r"}$=") + \
                       ("{:.3g}").format(q_num))
        plt.axvline(q_num,label=q_label,linestyle='--',linewidth=4)
    plt.xscale('log')
    plt.xlim([min(distance_limits),2])
    plt.ylim(0.5,max(plt.ylim()))
    PlotUtilities.lazyLabel(xlabel,"Count","",frameon=False,loc=loc)

def _gen_rupture_hist(to_bin,alpha=0.3,linewidth=0,**kwargs):
    """
    plots a generic rupture (or not, if there's no data)"

    Args:
        to_bin: what to histogram
        others: see plt.hist
    Returns: 
        nothing
    """
    if len(to_bin) == 0:
        return [],[],[]
    return plt.hist(to_bin,alpha=alpha,linewidth=linewidth,**kwargs)

def rupture_force_histogram(objs,**kwargs):
    """
    plots a rupture force histogram

    Args:
        obj: list of rupture objects
        others: see _gen_rupture_hist
    Returns: 
        nothing
    """
    ruptures,_ = Learning.get_rupture_in_pN_and_loading_in_pN_per_s(objs)
    return _gen_rupture_hist(ruptures,**kwargs)

def loading_rate_histogram(objs,**kwargs):
    """
    plots a loaifng rate histogram

    Args:
        obj: list of rupture objects
        others: see _gen_rupture_hist
    Returns: 
        nothing
    """
    _,loading_rates = Learning.get_rupture_in_pN_and_loading_in_pN_per_s(objs)
    return _gen_rupture_hist(loading_rates,**kwargs)

def rupture_plot(true,pred,fig,count_ticks=3,
                 scatter_kwargs=None,style_pred=None,
                 style_true=None,use_legend=True,count_limit=None,
                 color_pred=None,color_true=None,
                 bins_load=None,bins_rupture=None,
                 remove_ticks=True,lim_plot_load=None,lim_plot_force=None,
                 title="",distance_histogram=None,gs=None,
                 limit_percentile=True):
    if (distance_histogram is None):
        n_rows = 2
        n_cols = 2
        widths = [2,1]
        heights = [2,1]
        offset=0
    else:
        n_rows = 2
        n_cols = 3
        widths = [2,2,1]
        heights = [2,2]
        offset=1
    if (gs is None):
        gs = gridspec.GridSpec(n_rows,n_cols,width_ratios=widths,
                               height_ratios=heights)
    subplot_f = lambda x: plt.subplot(x)
    ruptures_true,loading_true = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(true)
    ruptures_pred,loading_pred = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(pred)
    if (color_pred is None):
        color_pred = color_pred_def
    if (color_true is None):
        color_true = color_true_def
    if (style_true is None):
        style_true = _style_true(color_true)
    if (style_pred is None):
        style_pred = _style_pred(color_pred)
    if (scatter_kwargs is None):
        scatter_kwargs = dict(style_true=dict(label="true",**style_true),
                              style_predicted=dict(label="predicted",
                                                   **style_pred))
    _lim_force,_bins_rupture,_lim_load,_bins_load = \
        Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                loading_true,loading_pred,
                                                limit=False)
    _lim_force_plot,_bins_rupture_plot,_lim_load_plot,_bins_load_plot = \
        Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                loading_true,loading_pred,
                                                limit=limit_percentile)
    if (bins_rupture is None):
        bins_rupture = _bins_rupture_plot
    if (bins_load is None):
        bins_load = _bins_load_plot
    if (lim_plot_load is None):
        lim_plot_load = _lim_load_plot
    if (lim_plot_force is None):
        lim_plot_force = _lim_force_plot
    if (distance_histogram is not None):
        ax_hist = plt.subplot(gs[:,0])
        histogram_event_distribution(**distance_histogram)
    ax0 = subplot_f(gs[0,offset])
    plot_true_and_predicted_ruptures(true,pred,**scatter_kwargs)
    PlotUtilities.xlabel("")
    plt.xlim(lim_plot_load)
    plt.ylim(lim_plot_force)
    PlotUtilities.title(title)
    if (remove_ticks):
        ax0.get_xaxis().set_ticklabels([])
    ax1 =subplot_f(gs[0,offset+1])
    hatch_true = true_hatch()
    true_style_histogram = _histogram_true_style(color_true=color_true,
                                                 label="true")
    pred_style_histogram = _histogram_predicted_style(color_pred=color_pred,
                                                     label="predicted")
    # for the rupture force, we dont add the label
    rupture_force_true_style = dict(**true_style_histogram)
    rupture_force_true_style['label'] = None
    rupture_force_pred_style = dict(**pred_style_histogram)
    rupture_force_pred_style['label'] = None
    rupture_force_histogram(pred,orientation='horizontal',bins=bins_rupture,
                            **rupture_force_pred_style)
    rupture_force_histogram(true,orientation='horizontal',bins=bins_rupture,
                            **rupture_force_true_style)
    PlotUtilities.lazyLabel("Count","","")
    ax = plt.gca()
    # push count to the top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top') 
    if (remove_ticks):
        ax1.get_yaxis().set_ticklabels([])
    if (count_limit is not None):
        plt.xlim(count_limit)
    plt.ylim(lim_plot_force)
    plt.xscale('log')
    ax4 = subplot_f(gs[1,offset])
    n_pred,_,_ = loading_rate_histogram(pred,orientation='vertical',
                                        bins=bins_load,
                                        **pred_style_histogram)
    n_true,_,_, = loading_rate_histogram(true,orientation='vertical',
                                         bins=bins_load,**true_style_histogram)
                                         
    if (count_limit is None and (len(n_pred) * len(n_true) > 0)):
        max_n = np.max([n_pred,n_true])
        count_limit = [0.5,max_n*10]
    else:
        count_limit = plt.ylim()
    PlotUtilities.lazyLabel("loading rate (pN/s)","Count","",frameon=False,
                            loc='upper left',useLegend=use_legend)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(lim_plot_load)
    plt.ylim(count_limit)
    ax3 = subplot_f(gs[1,offset+1])
    if (len(loading_pred) > 0):
        coeffs = Analysis.\
            bc_coeffs_load_force_2d(loading_true,loading_pred,bins_load,
                                    ruptures_true,ruptures_pred,bins_rupture)
        # just get the 2d (last one
        coeffs = [1-coeffs[-1]]
    else:
        coeffs = [0]
    labels_coeffs = [r"BCC"]
    # add in the relative distance metrics, if the are here
    if (distance_histogram is not None):
        _,_,cat_relative_median,cat_relative_q,q = \
            Offline.relative_and_absolute_median_and_q(**distance_histogram)
        coeffs.append(cat_relative_q)
        q_fmt = str(int(q))
        labels_coeffs.append(r"P$_{" + q_fmt + "}$")
    index = np.array([i for i in range(len(coeffs))])
    bar_width = 0.5
    rects1 = plt.bar(index, coeffs,alpha=0.3,color=color_pred)
    label_func = lambda i,r: "{:.3g}".format(r.get_height())
    y_func = lambda i,r: r.get_height()/2
    PlotUtilities.autolabel(rects1,label_func=label_func,y_func=y_func,
                            fontsize=PlotUtilities.g_font_legend,
                            fontweight='bold')
    plt.xticks(index, labels_coeffs,
               rotation=30,fontsize=PlotUtilities.g_font_label)
    PlotUtilities.ylabel("Metric")
    PlotUtilities.tickAxisFont()
    # push metric to the right
    ax = plt.gca()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right') 
    ax.tick_params(axis=u'x', which=u'both',length=0)
    plt.ylim([0,1])

def rupture_distribution_plot(learner,out_file_stem,distance_histogram=dict()):
    """
    plots *and saves* the distributions of ruptures in *all* validation folds
    (ie: the entire sample) 

    Args:
        learner: see plot_individual_learner
        out_file_stem: base to save
    Returns:
        nothing
    """
    name = learner.description.lower()
    x_values = learner.param_values()
    ruptures_valid_true,ruptures_valid_pred = \
        Learning.get_true_and_predicted_ruptures_per_param(learner)
    for i,(param,true,pred) in enumerate(zip(x_values,ruptures_valid_true,
                                             ruptures_valid_pred)):
        fig = PlotUtilities.figure(figsize=(16,8))
        rupture_plot(true,pred,fig,distance_histogram=distance_histogram,
                     limit_percentile=True)
        out_path = "{:s}{:s}{:d}.png".format(out_file_stem,name,i)
        PlotUtilities.savefig(fig,out_path)

def plot_individual_learner(cache_directory,learner,rupture_kwargs=dict()):
    """
    Plots the results for a single, individual learner

    Args:
        cache_directory: where to save the plots
        learner: learning_curve instance to use
    Returns:
        nothing
    """
    learner_name = learner.description
    out_file_stem = cache_directory + "{:s}".format(learner_name)
    # get the scoring objects by paramter by fold
    train_scores = learner._scores_by_params(train=True)
    valid_scores = learner._scores_by_params(train=False)
    x_values = learner.param_values()
    rupture_distribution_plot(learner,out_file_stem,
                              distance_histogram=rupture_kwargs)
    fig = PlotUtilities.figure()
    distance_distribution_plot(learner,to_true=True)
    PlotUtilities.savefig(fig,out_file_stem + "histogram_to_true.png")
    fig = PlotUtilities.figure()
    distance_distribution_plot(learner,to_true=False)
    PlotUtilities.savefig(fig,out_file_stem + "histogram_to_predicted.png")
    fig = PlotUtilities.figure()
    plot_num_events_off(x_values,train_scores,valid_scores)
    PlotUtilities.savefig(fig,out_file_stem + "n_off.png")
    fig = PlotUtilities.figure()
    cross_validation_distance_metric(x_values,train_scores,valid_scores,
                                     to_true=True)
    PlotUtilities.savefig(fig,out_file_stem + "dist.png")

def plot_no_event(x,y,interp,slice_fit,probability_distribution,stdev_masked,
                  sigma,epsilon):
    plt.subplot(3,1,1)
    plt.plot(x,y * 1e12)
    plt.plot(x,interp(x) * 1e12)
    PlotUtilities.lazyLabel("","force","")
    plt.subplot(3,1,2)
    plt.plot(x,stdev_masked,label="masked stdev")
    plt.axhline(epsilon+sigma,color='b')
    plt.axhline(epsilon,color='r',label="epsilon")
    plt.axhline(epsilon-sigma,color='b')
    PlotUtilities.lazyLabel("","stdev","")
    plt.subplot(3,1,3)
    plt.plot(x,probability_distribution)
    PlotUtilities.lazyLabel("Time","Probability","")
    plt.yscale('log')

def top_bars(x,y,slices,colors,ymin=None,ymax=None):
    """
    adds bars to show when each slice changes (by color) from ymin to ymax

    Args:
        x,y: plot units
        slices: array of slices: same size as colors
        colors: given to axvspan
        ymin/ymax: how thick the bar should be in y
    Returns:
        nothing
    """
    ylim = plt.ylim()
    if (ymin is None):
        ymin = 0.8
    if (ymax is None):
        ymax = 0.90
    for s,c in zip(slices,colors):
        x_sliced = x[s]
        plt.axvspan(xmin=x_sliced[0],xmax=x_sliced[-1],ymin=ymin,ymax=ymax,
                    color=c,alpha=0.3,linewidth=0)

def before_and_after(x,y,before_slice,after_slice,style=dict(),
                     color_before='k',color_after='r',label=None,
                     label_for_before=True):
    """
    plots x and y two before and after slices

    Args:
        x,y: the x and y to plot
        before_slice: what slice to plot before
        after_slice: what slice to plot after
        style; for each of them
        label: for one of them
    """
    if label_for_before:
        label_before = label
        label_after = None
    else:
        label_before = None
        label_after = label
    before_slice = slice(before_slice.start,before_slice.stop+1,1)
    tuples = [ [x,y,before_slice,color_before,style,label_before],
               [x,y,after_slice,color_after,style,label_after]]
    for x_tmp,y_tmp,slice_v,color_tmp,style_tmp,label in tuples:
        x_sliced = x_tmp[slice_v]
        plt.plot(x_sliced,y_tmp[slice_v],color=color_tmp,label=label,
                 **style_tmp)

def plot_fec(example,colors=_fec_event_colors,n_filter=1000,use_events=True):
    """
    plots the given fec (*not* split)

    Args:
        example: TimeSepForce to split 
        colors: if use_events, one color per event on the fec (we switch)
        n_filter: how many points to use while filtering
    Returns:
        None
    """
    fec_split = Analysis.zero_and_split_force_extension_curve(example)
    retract = fec_split.retract
    retract.Force -= np.median(retract.Force)
    retract_filtered = FEC_Util.GetFilteredForce(retract,n_filter)
    # get everything in terms of ploting variables
    x_plot = lambda x: x * 1e9
    y_plot = lambda y: y * 1e12
    sep = x_plot(retract.Separation)
    force = y_plot(retract.Force)
    sep_filtered = x_plot(retract_filtered.Separation)
    force_filtered = y_plot(retract_filtered.Force)
    if fec_split.has_events() and use_events:
        slices = fec_split.get_retract_event_slices()
        colors_before = colors
        colors_after = colors
        for i in range(len(slices)-1):
            before_kwargs = dict(before_slice=slices[i],after_slice=slices[i+1],
                                 color_before=colors_before[i],
                                 color_after=colors_after[i+1])
            before_and_after(x=sep,y=force,style=dict(alpha=0.3),
                             **before_kwargs)
            before_and_after(x=sep_filtered,y=force_filtered,
                             style=dict(alpha=1),**before_kwargs)
    else:
        style = dict(color=colors[0])
        plt.plot(sep,force,alpha=0.3,**style)
        plt.plot(sep_filtered,force_filtered,**style)
    return fec_split

def plot_arrows_above_events(event_idx,plot_x,plot_y,fudge_y,color='g',
                             marker='v',markersize=15,alpha=1,zorder=10,
                             **kwargs):
    """
    plots arrows at the given indices, signifying an event

    Args:
        event_idx: where to put the arrow as indices into plot_<x/y>
        fudge_y: array of offsets (ie: for moving the arrow up to prevent 
        obscuring data)
    
        others: see plt.plot()
    """
    kw = dict(facecolor=color,
              zorder=zorder,
              marker=marker,
              s=markersize,
              alpha=alpha,**kwargs)
    for start in event_idx:
        plt.scatter(plot_x[start],plot_y[start]+fudge_y,**kw)

def plot_format(time_sep_force):
    x_plot = time_sep_force.Time.copy()
    x_plot -= min(x_plot)
    y_plot = time_sep_force.Force.copy() * 1e12
    return x_plot,y_plot
