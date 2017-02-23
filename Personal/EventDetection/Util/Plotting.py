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
import Analysis

def plot_distribution(x,f,f_interp,f_cdf,thresh,num_bins=50):
    """
    plots a distribution
    
    Args:
        x: time / separation (whatever abscicca is)
        f,f_interp: the raw and interpolated f values
        f_cdf: the cdf we are using ; low probability means an event is likely
        thresh: the minimum probability
        num_bins: number of bins for showing the histogram (purely for plots)
    Returns:
        nothing, plots the distribution
    """
    idx_high = np.where(f_cdf >= thresh)
    idx_low = np.where(f_cdf <= thresh)
    events = f_cdf[idx_low]
    f_minus_mu = f - f_interp
    mu,std = Analysis.spline_residual_mean_and_stdev(f,f_interp)
    limits = [min(f_minus_mu),max(f_minus_mu)]
    num_bins = 50
    linspace_f_diff = np.linspace(*limits,num=num_bins)
    pdf_diff = norm.pdf(x=linspace_f_diff,loc=mu,scale=std)
    plt.subplot(3,1,1)
    plt.plot(x,f_interp,color='b',linewidth=3)
    plt.plot(x,f,color='k',alpha=0.3)
    PlotUtilities.lazyLabel("Time (au)","Force (au)","")
    plt.subplot(3,1,2)
    plt.hist(f_minus_mu,bins=num_bins,normed=True)
    plt.plot(linspace_f_diff,pdf_diff,linewidth=3,color='r')
    PlotUtilities.lazyLabel("Force Difference (au)","Probability (au)","")
    plt.subplot(3,1,3)
    plt.semilogy(x[idx_high],f_cdf[idx_high],alpha=0.3,color='k',
                 label="Non-Event")
    plt.semilogy(x[idx_low],f_cdf[idx_low],linestyle='None',marker='.',
                 color='r',label="Event")
    PlotUtilities.lazyLabel("Time (au)","Probability","",frameon=True,
                            loc="lower right")

def plot_surface_idx(surface_idx,n_smooth,Obj):
    smoothed = FEC_Util.GetFilteredForce(Obj,n_smooth)
    x,f = Obj.Time, Obj.Force
    x_smoothed, f_smoothed = smoothed.Time, smoothed.Force
    plt.plot(x,f,color='k',alpha=0.3)
    plt.plot(x_smoothed,f_smoothed,color='b')
    plt.plot(x[surface_idx],f[surface_idx],'ro')
                            
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
    event_idx_end,event_idx_start,event_idx = info.end,\
                                              info.start,\
                                              info.event_idx
    mask = info.mask
    interp_first_deriv = info.interp.derivative(1)(time)
    # get the interpolated derivative
    interpolated_force = info.interp(time)
    tau = ex.tau
    stdevs = info.local_stdev
    # plot everything
    style_events = dict(color='r',label="True events")
    n_plots = 3
    x = x_func(retract)
    min_x,max_x = min(x),max(x)
    ylabel = "Force [pN]"
    x_range = max_x - min_x
    fudge = x_range * 0.05
    x_limits = [min_x - fudge,max_x + fudge]
    force_plot = force * 1e12
    plt.subplot(n_plots,1,1)
    plt.plot(x,force_plot,color='k',alpha=0.3)
    plt.plot(x,interpolated_force,color='b',linewidth=2)
    highlight_events(event_slices,x,force_plot,**style_events)
    PlotUtilities.lazyLabel("",ylabel,"")
    plt.xlim(x_limits)
    plt.subplot(n_plots,1,2)
    # plot the autocorrelation time along the plot
    min_x_auto = min(x) * 1.1
    auto_correlation_x = [min_x_auto,min_x_auto+tau]
    plt.semilogy(x[info.slice_fit],info.cdf)
    highlight_events(event_slices,x,info.cdf,linewidth=5,**style_events)
    plt.axhline(thresh,label="threshold",linestyle='--',color='r')
    PlotUtilities.lazyLabel("","No-Event CDF ","")
    plt.xlim(x_limits)
    plt.subplot(n_plots,1,3)
    # XXX check mask has at least one...
    plt.plot(x,force_plot,'b-',color='k',alpha=0.3)
    for fwd,rev,event in zip(event_idx_end,event_idx_start,event_idx):
        plt.axvline(x[fwd],linestyle='--',color='r')
        plt.axvline(x[rev],color='g')
        plt.axvline(x[event],linewidth=3)
    plt.xlim(x_limits)
    PlotUtilities.lazyLabel(xlabel,ylabel,"")

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

    
    

