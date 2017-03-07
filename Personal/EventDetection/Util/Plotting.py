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
import Analysis,Learning

style_train = dict(color='r',marker='o',linestyle='--',label="Training") 
style_valid = dict(color='g',marker='v',linestyle='-',label="Validation")

                            
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
    plt.ylim([np.percentile(log_norm,0.5),max(log_norm)])
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
    interp_first_deriv = info.interp.derivative(1)(time)
    # get the interpolated derivative
    interpolated_force = info.interp(time)
    tau = ex.tau
    stdevs = info.local_stdev
    # plot everything
    style_events = dict(color='r',label="True events")
    x = x_func(retract)
    surface_index = ex.get_predicted_retract_surface_index()
    min_x,max_x = min(x),max(x)
    ylabel = "Force [pN]"
    x_range = max_x - min_x
    fudge = x_range * 0.05
    x_limits = [min_x - fudge,max_x + fudge]
    force_plot = force * 1e12
    interpolated_force_plot = interpolated_force*1e12
    # get the informaiton relevant to the CDF
    cdf = info.cdf
    boolean_mask = np.zeros_like(cdf)
    boolean_mask[mask] = 1
    masked_cdf = cdf.copy()
    masked_cdf *= boolean_mask
    n_rows = 4
    n_cols = 1
    lazy_kwargs = dict(frameon=True,loc="best")
    plt.subplot(n_rows,n_cols,1)
    plt.plot(x,force_plot,color='k',alpha=0.3,label="data")
    plt.plot(x,interpolated_force_plot,color='b',linewidth=2,label="2-spline")
    plt.axvline(x[surface_index],label="Predicted surface location")
    highlight_events(event_slices,x,force_plot,**style_events)
    PlotUtilities.lazyLabel("",ylabel,"",**lazy_kwargs)
    plt.xlim(x_limits)
    plt.subplot(n_rows,n_cols,2)
    # plot the autocorrelation time along the plot
    min_x_auto = min(x) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+tau]
    plt.semilogy(x,cdf,color='k',alpha=0.3,label="cdf")
    plt.semilogy(x,masked_cdf,color='b',linewidth=1,label="masked cdf")
    plt.axhline(thresh,color='k',linestyle='--',label="threshold")
    mask_boolean = np.zeros(x.size)
    mask_boolean[mask] = 1
    PlotUtilities.lazyLabel("","No-Event CDF ","",**lazy_kwargs)
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
    plt.subplot(n_rows,n_cols,4)
    mask_styles = [dict(linewidth=3,color='k',linestyle='-.',alpha=0.3),
                   dict(linewidth=1,color='r',linestyle='--',alpha=0.7)]
    tol = 1e-6
    for i,c in enumerate(info.condition_results):
        bool = np.zeros_like(x)
        bool[c] = 1
        style = mask_styles[i % len(mask_styles)]
        plt.semilogy(x,bool+tol,label="mask {:d}".format(i),**style)
    plt.plot(x,boolean_mask,color='k',linestyle='-',label="final mask")
    plt.xlim(x_limits)
    PlotUtilities.lazyLabel(xlabel,"mask (au)","",**lazy_kwargs)

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
    PlotUtilities.lazyLabel("","Force (pN)","",frameon=True)
    plt.subplot(2,1,2)
    plt.plot(x_plot(fit_x),y_plot(fit_y))
    plt.plot(x_plot(x_event),y_plot(y_event),color='r',alpha=0.3,label="event")
    plt.plot(x_plot(x_event),y_plot(pred),label="prediction")
    if (len(idx_above_predicted) > 0):
        plt.axvline(x_plot(x[idx_above_predicted[-1]]),label="surface pred")
    PlotUtilities.lazyLabel("Time","Force (pN)","",frameon=True,
                            loc="upper right")

def debug_plot_adhesion_info(probability_distribution,no_event_mask,min_idx,
                             no_event_boundaries,event_boundaries,threshold,
                             to_ret):
    """
    Used  insides of Detector.adhesion_mask to tell wtf is happening
    
    Args:
        all: internal Detector.adhesion_mask, see that function
    Returns:
        nothing
    """                             
    boolean_debug = np.zeros(probability_distribution.size)
    boolean_debug[no_event_mask] = 1
    plt.subplot(2,1,1)
    plt.plot(probability_distribution,color='g',linewidth=0.5)
    for i,e in enumerate(no_event_boundaries):
        label = "no event start" if i ==0 else None
        plt.axvline(e.start,linestyle='-',color='r',alpha=0.3,
                    label=label)
    for i,e in enumerate(event_boundaries):
        label_start = "event start" if i ==0 else None    
        label_end = "event end" if i ==0 else None  
        plt.axvline(e.start,linestyle='-.',color='b',alpha=0.5,
                    label=label_start)
        plt.axvline(e.stop,linestyle='-.',color='b',linewidth=4,alpha=0.5,
                    label=label_end)        
    plt.axvline(min_idx,linewidth=4,color='k',alpha=0.3,linestyle="--",
                label="minimum index used")
    plt.axhline(threshold,label="probability threshold")
    plt.yscale('log')
    PlotUtilities.lazyLabel("","probability","",loc="upper right",
                            frameon=True) 
    plt.subplot(2,1,2)
    plt.plot(to_ret)
    PlotUtilities.lazyLabel("index (au)","boolean no event mask","",
                            loc="upper right")
        
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
    to_pN = lambda x: x * 1e12
    rupture_forces_pN = [to_pN(obj.rupture_force) for obj in to_plot]
    loading_rate_pN_per_s = [to_pN(obj.loading_rate) for obj in to_plot]
    plt.semilogx(loading_rate_pN_per_s,rupture_forces_pN,**kwargs)

def plot_predicted_and_true_ruptures(true,predicted,title="",label_true="true",
                                     style_predicted=None):
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
    line_style = dict(linestyle="None")
    if (style_predicted is None):
        style_predicted  = dict(marker='x',color='k',label="predicted",
                                linewidth=2,**line_style)
    style_true = dict(marker='o',color='g',label=label_true,alpha=0.5,
                      linewidth=0,**line_style)
    _plot_rupture_objects(true,**style_true)
    _plot_rupture_objects(predicted,**style_predicted)
    PlotUtilities.lazyLabel("Loading Rate [pN/s]","Rupture Force [pN]",title,
                            frameon=True,legend_kwargs=dict(numpoints=1))



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
    fig = PlotUtilities.figure(figsize=(8,16))    
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
    plot_predicted_and_true_ruptures(rupture_true,rupture_predicted)

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
                            frameon=True)
    
def plot_num_events_off(x_values,train_scores,valid_scores):
    """
    Plots the number of 

    Args:
        cache_directory: where to save the plots
        learner: learning_curve instance to use
    Returns:
        nothing
    """
    x_train,train_dist,train_dist_std = \
        Learning.number_events_off_per_param(x_values,train_scores)
    x_valid,valid_dist,valid_dist_std = \
        Learning.number_events_off_per_param(x_values,train_scores)
    train_dist_plot,train_error_plot = train_dist,train_dist_std
    valid_dist_plot,valid_error_plot = valid_dist,valid_dist_std
    plt.errorbar(x=x_train,y=train_dist_plot,yerr=train_error_plot,
                 **style_train)
    plt.errorbar(x=x_valid,y=valid_dist_plot,yerr=valid_error_plot,
                 **style_valid)
    PlotUtilities.lazyLabel("Tuning parameter",
                            "Relative number of missing or incorrect events",
                            "")
    plt.xscale('log')    
    plt.yscale('log')    

def distance_distribution_plot(learner,**kwargs):
    """
    plots the distribution of distances to/from predicted events from/to
    actual events, dependning on kwargs
    
    Args:
        learner: the learner object to use
        kwargs: passed to event_distance_distribution (ie: to_true=T/F)
    """
    train_scores = learner._scores_by_params(train=True)
    valid_scores = learner._scores_by_params(train=False)
    name = learner.description.lower()
    x_values = learner.param_values()
    train_dist = Learning.event_distance_distribution(train_scores,**kwargs)
    valid_dist = Learning.event_distance_distribution(valid_scores,**kwargs)
    dist_plot = lambda x: [v * 1e9 for v in x]
    train_plot = dist_plot(train_dist)
    valid_plot = dist_plot(valid_dist)
    plt.boxplot(x=train_plot)
    plt.boxplot(x=valid_plot)
    plt.gca().set_yscale('log')
    PlotUtilities.lazyLabel("Tuning parameter","Distance Distribution (nm)",
                            "Event distributions for {:s}".format(name))
    

def plot_individual_learner(cache_directory,learner):
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
