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
import matplotlib.gridspec as gridspec

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
    original_cdf = info.probabilities[0]
    cdf = original_cdf
    boolean_mask = np.zeros_like(cdf)
    boolean_mask[mask] = 1
    masked_cdf = cdf.copy()
    masked_cdf *= boolean_mask
    n_rows = 3
    n_cols = 1
    lazy_kwargs = dict(frameon=True,loc="best")
    plt.subplot(n_rows,n_cols,1)
    plt.plot(x,force_plot,color='k',alpha=0.3,label="data")
    plt.plot(x,interpolated_force_plot,color='b',linewidth=2,label="2-spline")
    plt.axvline(x[surface_index],label="Surface\n(pred)")
    highlight_events(event_slices,x,force_plot,**style_events)
    PlotUtilities.lazyLabel("",ylabel,"",**lazy_kwargs)
    plt.xlim(x_limits)
    plt.subplot(n_rows,n_cols,2)
    # plot the autocorrelation time along the plot
    min_x_auto = min(x) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+tau]
    styles = [dict(color='k',linestyle='-',alpha=0.3),
              dict(color='g',linestyle='-.',alpha=0.7),
              dict(color='r',linestyle=':'),
              dict(color='m',linestyle='--')]
    for i,c in enumerate(info.probabilities):
        sty = styles[i % len(styles)]
        plt.semilogy(x,c,label="cdf{:d}".format(i),**sty)
    min_cdf = min([min(c) for c in info.probabilities])
    plt.semilogy(x,masked_cdf,color='b',linewidth=1,label="masked cdf")
    plt.axhline(thresh,color='k',linestyle='--',label="threshold")
    mask_boolean = np.zeros(x.size)
    mask_boolean[mask] = 1
    PlotUtilities.lazyLabel("","No-Event CDF ","",**lazy_kwargs)
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
    PlotUtilities.lazyLabel("","Force (pN)","",frameon=True)
    plt.subplot(2,1,2)
    plt.plot(x_plot(fit_x),y_plot(fit_y))
    plt.plot(x_plot(x_event),y_plot(y_event),color='r',alpha=0.3,label="event")
    plt.plot(x_plot(x_event),y_plot(pred),label="prediction")
    if (len(idx_above_predicted) > 0):
        plt.axvline(x_plot(x[idx_above_predicted[-1]]),label="surface pred")
    PlotUtilities.lazyLabel("Time","Force (pN)","",frameon=True,
                            loc="upper right")

def debug_plot_adhesion_info(time,force,force_fit,min_idx,derivative_gt_zero,
                             derivative_le_zero,to_ret):
    """
    Used inside of Detector.adhesion_mask to tell wtf is happening
    
    Args:
        all: internal Detector.adhesion_mask, see that function
    Returns:
        nothing
    """                             
    plt.subplot(2,1,1)
    plt.plot(time,force*1e12,color='k',alpha=0.3)
    plt.plot(time,force_fit*1e12)
    plt.axvline(time[min_idx])
    PlotUtilities.lazyLabel("","Force","",loc="upper right",
                            frameon=True)     
    plt.subplot(2,1,2)
    plt.plot(time,derivative_gt_zero)
    plt.plot(time,derivative_le_zero)
    plt.plot(time,to_ret,color='k',linestyle='--')
    PlotUtilities.lazyLabel("Time","mask","",loc="upper right",
                            frameon=True)   

        
def get_rupture_in_pN_and_loading_in_pN_per_s(objs):
    """
    Args:
        objs: see _plot_rupture_objecs
    Returns:
        tuple of <rupture force in pN, loading rate in pN>
    """
    to_pN = lambda x: x * 1e12
    rupture_forces_pN = np.array([to_pN(obj.rupture_force) for obj in objs])
    loading_rate_pN_per_s = np.array([to_pN(obj.loading_rate) for obj in objs])
    return rupture_forces_pN,loading_rate_pN_per_s

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
        get_rupture_in_pN_and_loading_in_pN_per_s(to_plot)
    plt.semilogx(loading_rate_pN_per_s,rupture_forces_pN,**kwargs)
    # XXX debugging...
    """
    from FitUtil.EnergyLandscapes.Lifetime_Dudko2008.Python.Code import \
        Dudko2008Lifetime
    forces = rupture_forces_pN * 1e-12
    loading_rates = loading_rate_pN_per_s *1e-12
    range_tau0 = [1e-3,1e3]
    range_x_tx = [10e-9,200e-9]
    kbT = 4.1e-21
    range_DeltaG_tx = [kbT,200*kbT]
    good_idx = np.where(loading_rates > 0)
    fit_dict = dict(ranges=[range_tau0,range_x_tx,range_DeltaG_tx],Ns=20)
    good_forces = forces[good_idx]
    fit = Dudko2008Lifetime.dudko_fit(good_forces,loading_rates[good_idx],
                                      fit_dict=fit_dict)
    space_forces = np.linspace(min(good_forces),max(good_forces))
    pred_rates = 1e12 * fit.predict(space_forces)
    plt.plot(pred_rates,space_forces*1e12,'r--',linewidth=3)
    print(pred_rates)
    plt.show()
    print("res!")
    print(fit.fit_result)
    exit(1)
    """

def plot_true_and_predicted_ruptures(true,predicted,title="",label_true="true",
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
    line_style = dict(linestyle="None")
    if (style_predicted is None):
        style_predicted = dict(color='k',label="predicted",
                                linewidth=2,**line_style)
    if (style_true is None):
        style_true = dict(color='g',label=label_true,alpha=0.5,**line_style)
    
    _plot_rupture_objects(true,marker='o',linewidth=0,linestyle="None",
                          **style_true)
    _plot_rupture_objects(predicted,marker='x',linewidth=3,linestyle="None",
                          **style_predicted)
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
                            "Event distributions for {:s}".format(name))
    
def _gen_rupture_hist(to_bin,alpha=0.3,linewidth=0,**kwargs):
    plt.hist(to_bin,alpha=alpha,linewidth=linewidth,**kwargs)

def rupture_force_histogram(objs,**kwargs):
    ruptures,_ = get_rupture_in_pN_and_loading_in_pN_per_s(objs)
    _gen_rupture_hist(ruptures,**kwargs)

def loading_rate_histogram(objs,**kwargs):
    _,loading_rates = get_rupture_in_pN_and_loading_in_pN_per_s(objs)
    _gen_rupture_hist(loading_rates,**kwargs)

def rupture_plot(true,pred,count_ticks=3,scatter_kwargs=None,style_pred=None,
                 style_true=None,
                 lim_load=None,lim_force=None,bins_load=None,bins_rupture=None,
                 remove_ticks=True):
    gs = gridspec.GridSpec(2, 2,
                           width_ratios=[4,1],
                           height_ratios=[4,1])
    ruptures_true,loading_true = \
        get_rupture_in_pN_and_loading_in_pN_per_s(true)
    ruptures_pred,loading_pred = \
        get_rupture_in_pN_and_loading_in_pN_per_s(pred)
    double_f = lambda f,*args: f([f(x) for x in args])
    if (style_true is None):
        style_true = dict(color='k',alpha=0.2)
    if (style_pred is None):
        style_pred = dict(color='g',alpha=0.7)
    if (scatter_kwargs is None):
        scatter_kwargs = dict(style_true=style_true,style_predicted=style_pred)
    if (lim_force is None):
        min_y = double_f(min,ruptures_pred,ruptures_true)
        max_y = double_f(max,ruptures_pred,ruptures_true)
        lim_force = [min_y/2,max_y*2]
    if (lim_load is None):
        safe = lambda x: [x[i] for i in np.where(np.array(x)>0)[0]]
        min_x = double_f(min,safe(loading_pred),safe(loading_true))
        max_x = double_f(max,safe(loading_pred),safe(loading_true))
        lim_load = [min_x*0.8,max_x*1.2]
    if (bins_rupture is None):
        bins_rupture= np.linspace(*lim_force,num=10)
    if (bins_load is None):
        min_y = max(min(lim_load),1e-2)
        logy = np.log10([min_y,max(lim_load)])
        bins_load = np.logspace(*logy,num=10)
    ax0 = plt.subplot(gs[0])
    plot_true_and_predicted_ruptures(true,pred,**scatter_kwargs)
    PlotUtilities.xlabel("")
    plt.xlim(lim_load)
    plt.ylim(lim_force)
    if (remove_ticks):
        ax0.get_xaxis().set_ticklabels([])
    ax1 = plt.subplot(gs[1])
    rupture_force_histogram(true,orientation='horizontal',bins=bins_rupture,
                            **style_true)
    rupture_force_histogram(pred,orientation='horizontal',bins=bins_rupture,
                            **style_pred)
    PlotUtilities.lazyLabel("Count","","")
    if (remove_ticks):
        ax1.get_yaxis().set_ticklabels([])
    plt.ylim(lim_force)
    ax4 = plt.subplot(gs[2])
    loading_rate_histogram(true,orientation='vertical',bins=bins_load,
                           **style_true)
    loading_rate_histogram(pred,orientation='vertical',bins=bins_load,
                           **style_pred)
    PlotUtilities.lazyLabel("loading rate [pN/s]","Count","")
    plt.xscale('log')
    plt.xlim(lim_load)
    ax3 = plt.subplot(gs[3])
    coeff_load = Analysis.\
                 bhattacharyya_probability_coefficient_1d(loading_true,
                                                         loading_pred,
                                                         bins_load)
    coeff_force = Analysis.\
                  bhattacharyya_probability_coefficient_1d(ruptures_true,
                                                          ruptures_pred,
                                                          bins_rupture)
    # do a 2-d coefficient
    tuple_true = [loading_true,ruptures_true]
    tuple_pred = [loading_pred,ruptures_pred]
    tuple_bins = [bins_load,bins_rupture]
    coeff_2d = Analysis.\
            bhattacharyya_probability_coefficient_dd(tuple_true,tuple_pred,
                                                     tuple_bins)
    coeffs = [coeff_load,coeff_force,coeff_2d]
    print(coeffs)
    labels_coeffs = [r"$\nu$",r"$F_r$",r"$\nu$,$F_r$"]
    index = np.array([i for i in range(len(coeffs))])
    bar_width = 0.5
    rects1 = plt.bar(index, coeffs,alpha=0.3,color='b')
    label_func = lambda i,r: "{:.2f}".format(r.get_height())
    y_func = lambda i,r: r.get_height()/2
    PlotUtilities.autolabel(rects1,label_func=label_func,y_func=y_func)
    plt.xticks(index + bar_width / 2, labels_coeffs,
               rotation=30,fontsize=PlotUtilities.g_font_legend)
    PlotUtilities.ylabel("BC value")
    # just empty :-(

def rupture_distribution_plot(learner,out_file_stem):
    """
    plots *and saves* the distributions of ruptures in *all* validation folds
    (ie: the entire sample) 

    Args:
        learner: see plot_individual_learner
        out_file_stem: base to save
    Returns:
        nothingS
    """
    train_scores = learner._scores_by_params(train=True)
    valid_scores = learner._scores_by_params(train=False)
    # get the validation ruptures (both truee and predicted)
    ruptures_valid_true = Learning.rupture_objects(valid_scores,get_true=True)
    ruptures_valid_pred = Learning.rupture_objects(valid_scores,get_true=False)
    name = learner.description.lower()
    x_values = learner.param_values()
    for i,(param,true,pred) in enumerate(zip(x_values,ruptures_valid_true,
                                             ruptures_valid_pred)):
        fig = PlotUtilities.figure()
        rupture_plot(true,pred)
        out_path = "{:s}{:s}{:d}.png".format(out_file_stem,name,i)
        PlotUtilities.savefig(fig,out_path)


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
    rupture_distribution_plot(learner,out_file_stem)
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

def debug_plot_force_value(x,f,interp_f,probability,probability_updated,
                           slice_to_use,bool_interp):
    """
    For debugging at the end of Detector.force_value_mask_function

    Args:
        see Detector.force_value_mask_function
    Returns: 
        nothing, makes a pretty plot.
    """
    force_plot = lambda x: x * 1e12
    plt.subplot(2,1,1)
    plt.plot(x, force_plot(f),alpha=0.3,color='k',label="raw")
    plt.plot(x, force_plot(interp_f),color='b',label="interpolated")
    PlotUtilities.lazyLabel("","Force","",loc='upper right')
    plt.subplot(2,1,2)
    plt.semilogy(x,bool_interp + 1e-2)
    plt.semilogy(x,probability_updated[slice_to_use],color='r',linestyle='--',
                 label="prob new")
    plt.semilogy(x,probability[slice_to_use],color='k',alpha=0.3,
                 label="prob orig")
    PlotUtilities.lazyLabel("Time","Prob/Mask","",loc='upper right')


def debug_plot_derivative(retract,slice_to_use,probability_updated,
                          boolean_ret,spline_probability_in_slice,
                          slice_updated,threshold,interp):
    """
    For debugging at the end of Detector.derivative_mask_function

    Args:
        see Detector.derivative_mask_function
    Returns: 
        nothing, makes a pretty plot.
    """
    time = retract.Time
    time_lim = [min(time),max(time)]
    x = retract.Time[slice_to_use]
    f = retract.Force * 1e12
    plt.subplot(2,1,1)
    plt.plot(time,f,label="force",color='k',alpha=0.3)
    plt.plot(x,f[slice_to_use])
    plt.plot(x,interp(x)*1e12,color='r')
    PlotUtilities.lazyLabel("","Force","",loc="upper right")
    plt.xlim(time_lim)
    plt.subplot(2,1,2)
    plt.plot(time,probability_updated,label="prob (updated)",
             linestyle='--')
    plt.plot(x,spline_probability_in_slice,color='k',alpha=0.3,
             label="prob (original)")
    plt.xlim(time_lim)
    plt.plot(time,boolean_ret+min(probability_updated),linewidth=1,alpha=0.3,
             label="mask")
    plt.axhline(threshold,label="t={:.3g}".format(threshold))
    mask_ends = np.zeros(time.size)
    mask_ends[slice_updated] = 1
    plt.plot(time,mask_ends,label="end mask")
    plt.yscale('log')
    plt.ylim([min(probability_updated)/5,2])
    PlotUtilities.lazyLabel("Time","Probability","",loc="upper right")
