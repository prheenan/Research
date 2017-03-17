# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,random,multiprocessing,re
from timeit import default_timer as timer
from GeneralUtil.python import GenUtilities,PlotUtilities

class timing_info:
    def __init__(self,learner_trials):
        n = len(learner_trials.list_of_time_trials)
        arr = lambda : np.zeros(n)
        self.nums,self.means,self.stdevs = [],[],[]
        self.velocities,self.pts_per,self.pts_std = arr(),arr(),arr()
        trials = learner_trials.list_of_time_trials
        # sort the trials by number of points
        pts = [t.average_number_of_points_per_curve() for t in trials]
        sort_idx = np.argsort(pts)
        trials = [trials[i] for i in sort_idx]
        for i,trial in enumerate(trials):
            num_curves = trial.num_curves
            mean = trial.mean_time_for_fixed_number()
            std = trial.std_time_for_fixed_number()
            self.nums.append(num_curves)
            self.means.append(mean)
            self.stdevs.append(std)
            self.velocities[i] = learner_trials.loading_rates[i]
            self.pts_per[i] = trial.average_number_of_points_per_curve()
            self.pts_std[i] = trial.stdev_number_of_points_per_curve()
    @property
    def size(self):
        return self.nums.size                           
                    


def get_loading_rate_slopes_and_errors(inf):
    """
    given a timing_info object, gets the slopes of runtime versus number of
    curves for each loading rate

    Args:
        inf: timing info object
    Returns:
        tuple of <slopes, slope errors>
    """
    coeffs = []
    for num,mean in zip(inf.nums,inf.means):
        coeffs.append(GenUtilities.GenFit(x=num,y=mean))
    params = [c[0][0] for c in coeffs]
    params_std = [c[1][0] for c in coeffs]
    return params,params_std

def get_linear_runtime(x,times,fudge=1):
    """
    given x and times (assumed d), gets the coefficients for a linear fit
    to times vs x

    Args:
        x: the x values
        times: time y values used
        fudge: when returning the prediction, how much to 'pad' in +/- x
    Returns:
        tuple of <coefficents for linear fit,their errors x predictions,
        y predictions>
    """
    sort_idx = np.argsort(x)
    x_pred = np.array(x)[sort_idx]
    y_pred = np.array(times)[sort_idx]
    params,params_std,_ = GenUtilities.GenFit(x=x_pred,y=y_pred)
    log_bounds = np.log10([min(x_pred)/fudge,max(x_pred)*fudge])
    x_pred_plot = np.logspace(*log_bounds,base=10,endpoint=True)
    y_pred_plot = np.polyval(params,x=x_pred_plot)
    return params,params_std,x_pred_plot,y_pred_plot


def plot_learner_prediction_time_comparison(learners,color='b'):
    """
    plots the asympotic slope of the learners

    Args:
        learners: list of learners
        color: passed to histogram
    Returns:
        nothing, makes a pretty plot
    """
    time_per_point,time_per_point_error = [],[]
    for l in learners:
        inf = timing_info(l)
        params,_ = get_loading_rate_slopes_and_errors(inf)
        x,_ = _timing_plot_pts_and_pts_error(inf)
        coeffs,coeffs_err,_,_ = get_linear_runtime(x,params)
        time_per_point.append(coeffs[0])
        time_per_point_error.append(coeffs_err[0])
    plot_y = 1/np.array(time_per_point)
    # d(1/f) = (1/f^2) * df
    plot_y_error = plot_y**2 * (np.array(time_per_point_error))
    N = len(time_per_point)
    labels = [l.learner.description.lower() for l in learners]
    ind = np.arange(N)  # the x locations for the groups
    width = 0.4       # the width of the bars
    ax = plt.gca()
    rects = ax.bar(ind , plot_y, width,color=color,alpha=0.4,linewidth=0,
                   yerr=plot_y_error,log=True,
                   error_kw=dict(ecolor='k',linewidth=2,capsize=15))
    # add some text for labels, title and axes ticks
    ax.set_xticks(ind)
    ax.set_xticklabels(labels)
    xlim = plt.xlim()
    fudge = width/2
    xlim = [min(xlim)-fudge,max(xlim) + fudge]
    plt.xlim(xlim)
    formatted_with_errors = [pretty_exp_with_error(r.get_height(),e) \
                             for r,e in zip(rects,plot_y_error)]
    # add labels with the errors for each
    label_func = lambda i,r : formatted_with_errors[i]
    kwargs = dict(fontsize=PlotUtilities.g_font_legend)
    PlotUtilities.autolabel(rects,label_func,fontdict=kwargs)
    PlotUtilities.lazyLabel("Event finding method",
                            "Points classified per second","")
    
        
def plot_learner_slope_versus_loading_rate(learner_trials,style_data=None,
                                           style_pred=None,min_pred=1e-3):
    """
    Makes a plot of the (slope of runtime versus number of curves) versus
    loading rate

    Args:
        learner_trials: a single learner object
        style_<data/pred>: style for the data points and predictions
    Returns:
        nothing, makes a pretty plot
    """
    if (style_data is None):
        style_data = dict(color='r',linestyle="None",marker='o')
    if (style_pred is None):
        style_pred = dict(color='b',linestyle="--")
    inf = timing_info(learner_trials)
    params,params_std = get_loading_rate_slopes_and_errors(inf)
    # the slope is the time per force extension curve (less an offset; get that
    # per loading rate
    velocities = inf.velocities
    x,xerr = _timing_plot_pts_and_pts_error(inf)
    coeffs,coeffs_err,x_pred,y_pred = get_linear_runtime(x,params)
    # fit a linear model to the runtime
    fudge = 1.75
    plt.errorbar(x=x,xerr=xerr,y=params,yerr=params_std,**style_data)
    slope = coeffs[0]
    lower_label = r"{:.2f}$ c_0 N \leq $".format(fudge)
    upper_label = r"$\leq \frac{c_0}{" + "{:.2f}".format(fudge) + "} N$"
    label_timing = learner_trials.learner.description.lower()
    style_timing = style_pred
    idx_good  = np.where(y_pred > min_pred)
    x_pred_plot = x_pred[idx_good]
    y_pred_plot = y_pred[idx_good]
    # only plot where the prediction is reasonable
    plt.plot(x_pred_plot,y_pred_plot/fudge,**style_timing)
    plt.plot(x_pred_plot,y_pred_plot*fudge,**style_timing)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    # plot something just for the legend entry
    plt.plot([],[],label=label_timing,marker=style_data['marker'],
             **style_timing)
    PlotUtilities.lazyLabel("N (points per curve)",
                            "Runtime per curve (s)",
                            "Runtime per curve, T(N), is $\Theta(\mathrm{N})$")



def _timing_plot_pts_and_pts_error(inf,factor=1):
    """
    gets the average number of points per curve in a given loading rate
    and the error (1 standard deviaiton)

    Args:
        inf: the timing_info object
        factor: what to divide by (defaults to thousand of points)
    Returns:
        tuple of <mean number of points per curve, stdev of points per curve>
    """
    n = inf.pts_per.size
    pts,pts_err = np.zeros(shape=n),np.zeros(shape=n)
    for i,(pts_tmp,xerr_tmp) in enumerate(zip(inf.pts_per,inf.pts_std)):
        pts[i] = pts_tmp/factor
        pts_err[i] = xerr_tmp/factor
    return pts,pts_err

def get_sigfig_sign_and_exponent(number,format_str="{:3.0e}"):
    """
    gets the significant figure(s), sign, and exponent of a number

    Args:
        number: the number we want
        format_str: how it should be formatted (limiting number of sig figs)
    Returns:
        tuple of <sig figs,sign,exponent>
    """
    scientific = format_str.format(number)
    pattern = r"""
               (\d+[\.]*\d*) # number.numbers
               e          # literal e
              ([+-])0*(\d+)     # either plus or minus, then exponent
              """
    sig = re.match(pattern,scientific,re.VERBOSE)
    return sig.groups()

def pretty_exp_with_error(number,error,error_fmt="{:.1f}",**kwargs):
    """
    retrns

    Args:
        number: the number we want
        error_str: how it should be formatted (limiting number of sig figs)
        **kwargs: passed to get_sigfig_sign_and_exponent
    Returns:
        pretty-printed (latex) of number +/- error, like: '(a +/- b) * 10^(c)'
    """
    sigfig,sign,exponent = get_sigfig_sign_and_exponent(number,**kwargs)
    # get the error in terms of the exponent of the number
    exponent_num = float(exponent) * -1 if sign == "-" else float(exponent)
    error_rel = error/(10**(exponent_num))
    string_number_and_error = sigfig + r"\pm" + error_fmt.format(error_rel)
    # add parenths
    string_number_and_error = "(" + string_number_and_error + ")"
    return _pretty_format_exp(string_number_and_error,sign,exponent)


def _pretty_format_exp(sig_fig,sign,exponent):
    """
    pretty prints the number sig_fig as <number * 10^(exponent)>
    
    Args:
        sig_fig: number to print
        sign: literal +/-
        exponent: what to put in 10^{right here}
    Returns:
        formatted string
    """
    sign_str = "" if sign == "+" else "-"
    to_ret = r"$" +  sig_fig + r"\cdot 10^{" + sign_str + exponent + r"}$" 
    return to_ret

def pretty_exp(number,**kwargs):
    """
    takes a number and returns its pretty-printed exponent format

    Args:
        number: see pretty_exp_with_error
        **kwargs: passed to get_sigfig_sign_and_exponent
    Returns:
        see pretty_exp_with_error
    """
    args = get_sigfig_sign_and_exponent(number,**kwargs)
    return _pretty_format_exp(*args)

def plot_learner_versus_loading_rate_and_number(learner_trials):    
    """
    makes a plot of the runtimes versus number of force extension curves
    for each loading rate used.
    
    Args:
        learner_trials: a single learner object
    Returns:
        nothing, makes a pretty plot
    """
    styles = [dict(color='r',marker='x',linestyle='--'),
              dict(color='b',marker='o',linestyle='-'),
              dict(color='k',marker='v',linestyle='-.'),
              dict(color='g',marker='s',linestyle='-',linewidth=3),
              dict(color='m',marker='*',linestyle='-.',linewidth=3),
              dict(color='k',marker='8',linestyle='-',linewidth=2),
              dict(color='r',marker='p',linestyle='-.',linewidth=2),
              dict(color='b',marker='d',linestyle='--'),
              dict(color='k',marker='<',linestyle='-'),
              dict(color='g',marker='+',linestyle='-.')]
    inf = timing_info(learner_trials)
    pts,xerr = _timing_plot_pts_and_pts_error(inf,factor=1)
    for i,(num,mean,yerr,vel) in enumerate(zip(inf.nums,inf.means,inf.stdevs,
                                               inf.velocities)):
        style = styles[i % len(styles)]
        velocity_label = r"v={:4d}nm/s".format(int(vel))
        number_pretty = r"N$\approx$"  + pretty_exp(pts[i])
        label = "{:s}\n({:s})".format(velocity_label,number_pretty)
        plt.errorbar(x=num,y=mean,yerr=yerr,label=number_pretty,alpha=0.7,
                     **style)
    title = "Runtime verus loading rate and number of curves\n" + \
           "(N, points/curve, in parenthesis) "
    PlotUtilities.lazyLabel("Number of Force-Extension Curves","Time",title)

