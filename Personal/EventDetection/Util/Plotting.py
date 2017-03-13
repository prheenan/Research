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
    styles = [dict(color='r',linestyle=':'),
              dict(color='k',linestyle='-.',alpha=0.3),
              dict(color='g',linestyle='-',alpha=0.7),
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
        plt.axvline(x_plot(x[idx_above_predicted[-1]]))
    PlotUtilities.lazyLabel("Time","Force (pN)","",frameon=True,
                            loc="upper left")

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
    plt.plot(time,derivative_gt_zero,label="ge")
    plt.plot(time,derivative_le_zero,label="le")
    plt.plot(time,to_ret,color='k',linestyle='--')
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

def plot_true_and_predicted_ruptures(true,predicted,title="",
                                     use_legend=True,loc='upper left',
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
        style_true = dict(color='g',label="true",alpha=0.5,**line_style)
    _plot_rupture_objects(true,marker='o',linewidth=0,linestyle="None",
                          **style_true)
    _plot_rupture_objects(predicted,marker='x',linewidth=3,linestyle="None",
                          **style_predicted)
    PlotUtilities.lazyLabel("Loading Rate [pN/s]","Rupture Force [pN]",title,
                            frameon=True,legend_kwargs=dict(numpoints=1),
                            useLegend=use_legend,loc=loc)



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
    
def plot_num_events_off(x_values,train_scores,valid_scores,ylim=None):
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
    if ylim is None:
        ylim = [1e-2,1]
    plt.ylim(ylim)
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

def distance_f_score_plot(distance_scores,bins=None,xlim_plot=None):
    """
    plots the distnce f scores 

    Args: 
        distance_scores: to plot 
        bins: how to bin the data. default to a huge logspace bining
    Returnss:
        noting
    """
    xlim = [np.min(distance_scores),np.max(distance_scores)]
    if (xlim_plot is None):
        xlim_plot = xlim
    if (bins is None):
        xlim_log = np.log10(xlim)
        bins = np.logspace(*xlim_log,num=10)
    plt.hist(distance_scores,bins=bins,log=True)
    plt.xscale('log')
    plt.xlim(xlim)
    PlotUtilities.lazyLabel("Count","Distance F-score","")

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
        return
    plt.hist(to_bin,alpha=alpha,linewidth=linewidth,**kwargs)

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
    _gen_rupture_hist(ruptures,**kwargs)

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
    _gen_rupture_hist(loading_rates,**kwargs)

def rupture_plot(true,pred,count_ticks=3,scatter_kwargs=None,style_pred=None,
                 style_true=None,use_legend=True,count_limit=None,
                 lim_load=None,lim_force=None,bins_load=None,bins_rupture=None,
                 remove_ticks=True,lim_plot_load=None,lim_plot_force=None):
    gs = gridspec.GridSpec(2, 2,
                           width_ratios=[2,1],
                           height_ratios=[2,1])
    ruptures_true,loading_true = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(true)
    ruptures_pred,loading_pred = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(pred)
    if (style_true is None):
        style_true = dict(color='k',alpha=0.5)
    if (style_pred is None):
        style_pred = dict(color='g',alpha=0.3)
    if (scatter_kwargs is None):
        scatter_kwargs = dict(style_true=dict(label="true",**style_true),
                              style_predicted=dict(label="predicted",
                                                   **style_pred))
    _lim_force,_bins_rupture,_lim_load,_bins_load = \
        Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                loading_true,loading_pred)
    if (lim_force is None):
        lim_force = _lim_force
    if (lim_load is None):
        bins_rupture = _bins_rupture
    if (bins_rupture is None):
        bins_rupture = _bins_rupture
    if (bins_load is None):
        bins_load = _bins_load
    if (lim_plot_load is None):
        lim_plot_load = lim_load
    if (lim_plot_force is None):
        lim_plot_force = lim_force
    ax0 = plt.subplot(gs[0])
    plot_true_and_predicted_ruptures(true,pred,use_legend=use_legend,
                                     **scatter_kwargs)
    PlotUtilities.xlabel("")
    plt.xlim(lim_plot_load)
    plt.ylim(lim_plot_force)
    if (remove_ticks):
        ax0.get_xaxis().set_ticklabels([])
    ax1 = plt.subplot(gs[1])
    hatch_true = "//"
    true_style_histogram = dict(hatch=hatch_true,zorder=10,**style_true)
    pred_style_histogam = dict(zorder=1,**style_pred)
    rupture_force_histogram(pred,orientation='horizontal',bins=bins_rupture,
                            **pred_style_histogam)
    rupture_force_histogram(true,orientation='horizontal',bins=bins_rupture,
                            **true_style_histogram)
    PlotUtilities.lazyLabel("Count","","")
    if (remove_ticks):
        ax1.get_yaxis().set_ticklabels([])
    if (count_limit is not None):
        plt.xlim(count_limit)
    plt.ylim(lim_plot_force)
    plt.xscale('log')
    ax4 = plt.subplot(gs[2])
    loading_rate_histogram(pred,orientation='vertical',bins=bins_load,
                          label="predicted", **pred_style_histogam)
    loading_rate_histogram(true,orientation='vertical',bins=bins_load,
                           label="true",**true_style_histogram)
    PlotUtilities.lazyLabel("loading rate [pN/s]","Count","",frameon=True,
                            loc='upper left',useLegend=use_legend)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(lim_plot_load)
    if (count_limit is not None):
        plt.ylim(count_limit)
    ax3 = plt.subplot(gs[3])
    if (len(loading_pred) > 0):
        coeffs = Analysis.\
            bc_coeffs_load_force_2d(loading_true,loading_pred,bins_load,
                                    ruptures_true,ruptures_pred,bins_rupture)
    else:
        coeffs = [0,0,0]
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
    plt.ylim([0,1])
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
    name = learner.description.lower()
    x_values = learner.param_values()
    ruptures_valid_true,ruptures_valid_pred = \
        Learning.get_true_and_predicted_ruptures_per_param(learner)
    for i,(param,true,pred) in enumerate(zip(x_values,ruptures_valid_true,
                                             ruptures_valid_pred)):
        fig = PlotUtilities.figure(figsize=(12,12))
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
    plt.semilogy(x,bool_interp + 1e-2,label="mask")
    plt.semilogy(x,probability_updated[slice_to_use],color='r',linestyle='--',
                 label="prob new")
    plt.semilogy(x,probability[slice_to_use],color='k',alpha=0.3,
                 label="prob orig")
    PlotUtilities.lazyLabel("Time","Prob/Mask","",loc='upper right')

def debug_plot_derivative_ratio(time,slice_to_use,
                                 ratio,interp_sliced,force_sliced,
                                 interp_slice_deriv,
                                 boolean_ret,probability_updated,
                                 absolute_min_idx,ratio_min_threshold):
    x_sliced = time[slice_to_use]
    xlim = [min(x_sliced),max(x_sliced)]
    where_possible = np.where(ratio < ratio_min_threshold)
    plot_interp_deriv = interp_slice_deriv/max(interp_slice_deriv)
    plt.subplot(3,1,1)
    plt.plot(x_sliced,interp_sliced*1e12,linewidth=3,label="interp")
    plt.plot(x_sliced,force_sliced*1e12,color='k',alpha=0.3,label="force")
    PlotUtilities.lazyLabel("","Force","")
    plt.xlim(xlim)
    plt.subplot(3,1,2)
    plt.plot(x_sliced,ratio,color='k',alpha=0.3,label="ratio") 
    plt.plot(x_sliced[where_possible],ratio[where_possible],color='b')
    PlotUtilities.lazyLabel("","ratio df/epsilon","")
    plt.axvline(time[absolute_min_idx])
    plt.xlim(xlim)
    plt.subplot(3,1,3)
    plt.plot(time,probability_updated,label="prob")
    plt.plot(time,boolean_ret + min(probability_updated),label="mask")
    plt.yscale('log')
    plt.xlim(xlim)
    PlotUtilities.lazyLabel("time","prob,mask","")

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

def debug_plot_derivs(approach_time,approach_force,
                      approach_interp_sliced,x_sliced,
                      force_sliced,interp_sliced,
                      approach_interp_deriv,interp_slice_deriv,
                      min_deriv):
    ylim = [min(force_sliced),max(force_sliced)]
    min_v = min([min(approach_interp_deriv),min(interp_slice_deriv)])
    max_v = max([max(approach_interp_deriv),max(interp_slice_deriv)])
    ylim_deriv = [min_v,max_v]
    plt.subplot(2,2,1)
    plt.plot(approach_time,approach_force,alpha=0.3)
    plt.plot(approach_time,approach_interp_sliced)
    plt.ylim(ylim)
    PlotUtilities.lazyLabel("time","Force","")
    plt.subplot(2,2,3)
    plt.plot(x_sliced,force_sliced,alpha=0.3)
    plt.plot(x_sliced,interp_sliced)
    plt.ylim(ylim)
    PlotUtilities.lazyLabel("time","Force","")
    plt.subplot(2,2,2)
    plt.plot(approach_time,approach_interp_deriv)
    plt.ylim(ylim_deriv)
    PlotUtilities.lazyLabel("time","Deriv","")
    plt.subplot(2,2,4)
    plt.plot(x_sliced, interp_slice_deriv)
    plt.axhline(min_deriv,label="Minimum of approach")
    plt.ylim(ylim_deriv)
    PlotUtilities.lazyLabel("time","Deriv","")
    

def before_and_after(x,y,before_slice,after_slice,style,label=None):
    """
    plots x and y two before and after slices

    Args:
        x,y: the x and y to plot
        before_slice: what slice to plot before
        after_slice: what slice to plot after
        style; for each of them
        label: for one of them
    """
    color_before = 'b'
    color_after = 'r'
    tuples = [ [x,y,before_slice,color_before,style,label],
               [x,y,after_slice,color_after,style,None]]
    for x_tmp,y_tmp,slice_v,color_tmp,style_tmp,label in tuples:
        x_sliced = x_tmp[slice_v]
        plt.plot(x_sliced,y_tmp[slice_v],color=color_tmp,label=label,
                 **style_tmp)

def debug_plot_derivative(retract,slice_to_use,probability_updated,
                          boolean_ret,probability_original,
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
    plt.plot(time,probability_original,color='k',alpha=0.3,
             label="prob (original)")
    plt.xlim(time_lim)
    plt.plot(time,boolean_ret+min(probability_original),linewidth=1,alpha=0.3,
             label="mask")
    plt.axhline(threshold,label="t={:.3g}".format(threshold))
    mask_ends = np.zeros(time.size)
    mask_ends[slice_updated] = 1
    plt.plot(time,mask_ends,label="end mask")
    plt.yscale('log')
    plt.ylim([min(probability_updated)/5,2])
    PlotUtilities.lazyLabel("Time","Probability","",loc="upper right")
