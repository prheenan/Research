# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,random,multiprocessing,re
from timeit import default_timer as timer

sys.path.append("../../../../")
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import Learning,InputOutput

class time_trials:
    def __init__(self,times,num_curves,fec_num_points):
        self.times = times
        self.num_curves = num_curves
        self.fec_num_points = fec_num_points
    def mean_time_for_fixed_number(self):
        """
        Returns:
            mean times to complete each of self.num_curves
        """
        return np.array([np.mean(x) for x in self.times])
    def std_time_for_fixed_number(self):
        """
        Returns:
            standard deviation of times to complete each of self.num_curves
        """    
        return np.array([np.std(x) for x in self.times])
    def total_number_of_points_per_curves(self):
        """
        Returns:
            total number of points predicted for each of self.num_curves
        """    
        return np.array([sum(n) for n in self.fec_num_points])

    def stdev_number_of_points_per_curves(self):
        """
        Returns:
            total number of points predicted for each of self.num_curves
        """    
        return [sum(n) for n in self.fec_num_points]

    def average_number_of_points_per_curve(self):
        """
        Returns:
            average number of points per curve over all curves
        """        
        return np.mean(np.concatenate(self.fec_num_points))
    def stdev_number_of_points_per_curve(self):
        """
        Returns:
            stdev number of points per curve over all curves
        """        
        return np.std(np.concatenate(self.fec_num_points))

    

class time_trials_by_loading_rate:
    def __init__(self,learner,list_of_time_trials,loading_rates):
        self.learner = learner
        # sort the time trials, sorting increasing by number of points
        idx_sort = np.argsort([t.average_number_of_points_per_curve
                               for t in list_of_time_trials])
        self.list_of_time_trials = [list_of_time_trials[i] for i in idx_sort]
        self.loading_rates = [loading_rates[i] for i in idx_sort]
    def max_time_trial(self):
        """
        Returns:
            the maximum time across all trials. useful for plotting
        """
        return max([max(np.concatenate(l.times))
                    for l in self.list_of_time_trials])
    def min_time_trial(self):
        """
        Returns:
            the minmumx time across all trials. useful for plotting
        """    
        return min([min(np.concatenate(l.times))
                    for l in self.list_of_time_trials])        


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
                    
def time_example_running(func,d):
    start = timer()
    func(d)
    elapsed_time =  timer() - start
    return elapsed_time
           
def time_example_multiproc(args):
    return time_example_running(*args)
             
def time_single(func,data,pool_size=None):
    """
    time a a single predicton of a set of data, per:
stackoverflow.com/questions/7370801/measure-time-elapsed-in-python/25823885#2582388

    Args:
        func: to call, just takes in the FEC
        data: what to use
    Returns:
        time, in seconds, that it takes to run. 
    """
    if (pool_size is None):
        pool_size = multiprocessing.cpu_count() - 1
    pool = multiprocessing.Pool(pool_size)   
    args = [ (func,d) for d in data]
    time_increments =  pool.map(time_example_multiproc,args)
    pool.close()    
    pool.join()        
    return sum(time_increments)
    
def get_single_trial_times(trials_per_curve_set,learner,data_tmp):
    times =[time_single(learner.func_to_call,data_tmp) 
            for t in range(trials_per_curve_set)]
    return times
    
def get_all_times(learner,data,list_of_curve_numbers,velocity,
                  cache_directory,force_trials,timing_threshold=None,
                  trials_per_curve_set=5):
    """
    gets the time_trials for a single loading rate (or whatever)

    Args:
         see cache_all_learners, except a single category's data
    Returns:
        a single time_trials object
    """
    if (timing_threshold is None):
        timing_threshold = 300
    times_all_trials = []
    sizes_all_trials = []
    num_curves_all_trials = []
    for i,l in enumerate(list_of_curve_numbers):
        # dont do trials that we cant actually time
        assert l <= len(data) , "Only {:d}, not {:d}, curves loaded".\
            format(l,len(data))
        # determine the data set we will use for this one
        data_tmp = data[:l]
        sizes = [d.Force.size for d in data_tmp]
        cache_id = "{:s}{:s}_v={:.1f}_l={:d}.pkl".\
                format(cache_directory,learner.description,velocity,l)
        times = CheckpointUtilities.\
            getCheckpoint(cache_id,get_single_trial_times,force_trials,
                          trials_per_curve_set,learner,data_tmp)
        times_all_trials.append(times)
        sizes_all_trials.append(sizes)
        num_curves_all_trials.append(len(data_tmp))
        # give up if it is taking too long, per curve
        average_time_per_curve = np.mean(times)
        if (average_time_per_curve > timing_threshold):
            break
    return time_trials(times_all_trials,num_curves_all_trials,sizes_all_trials)

def single_learner(learner,curve_numbers,categories,**kwargs):
    """
    gets the time_trials_by_loading_rate for a single object

    Args:
         see cache_all_learners, except a single learner
    Returns:
        a single time_trials_by_loading_rate
    """
    trials = []
    # get the time trials for each category (loading rate)
    loading_rates = [c.velocity_nm_s for c in categories]    
    for c,velocity in zip(categories,loading_rates):
        data = c.data
        trials.append(get_all_times(learner,data,curve_numbers,
                                    velocity=velocity,**kwargs))
    return time_trials_by_loading_rate(learner,trials,loading_rates)

def cache_all_learners(learners,categories,curve_numbers,cache_directory,
                       force=True,**kwargs):
    """
    caches and returns the timing results for all the learners

    Args:
        learners: list ofLearning.learnin_curve object to use
        categories: list of Learninrg.ForceExtensionCtegories to use
        cache_directory: where the cache lives
        curve_numbers: how many curves to use (from each category;
        they are all used separately)
        **kwargs: passed to curve_numbers
    Returns:
        list of time_trials_by_loading_rate objects
    """
    times = []
    # read in all the data
    categories = InputOutput.read_categories(categories,force,cache_directory,
                                             limit=max(curve_numbers))
    # get all the trials for all the learners        
    for l in learners:
        learner_file = (cache_directory + "_l_" + l.description)
        t = CheckpointUtilities.\
            getCheckpoint(learner_file,single_learner,
                          force,l,curve_numbers,categories,
                          cache_directory=cache_directory,
                          force_trials=force,**kwargs)
        times.append(t)
    return times
        
        
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    learners = Learning.get_learners()
    positives_directory = InputOutput.get_positives_directory()
    positive_categories = InputOutput.\
        get_categories(positives_directory=positives_directory,
                       use_simulated=True)
    curve_numbers = [1,2,5,10,30,50,100,150,200]
    cache_dir = "../_1ReadDataToCache/cache/"
    force = False
    times = CheckpointUtilities.getCheckpoint(cache_dir + "all.pkl",
                                              cache_all_learners,force,
                                              learners,positive_categories,
                                              curve_numbers,
                                              cache_dir,force)
    # sort the times by their loading rates
    max_time = max([l.max_time_trial() for l in times])
    min_time = min([l.min_time_trial() for l in times])
    # plot the Theta(n) coefficient for each
    fig = PlotUtilities.figure()
    plot_learner_prediction_time_comparison(times)
    PlotUtilities.legend(loc="lower right",frameon=True)
    PlotUtilities.savefig(fig,"compare.png")
    for learner_trials in times:
        base_name = learner_trials.learner.description
        # plot the timing veruses loading rate and number of points 
        fig = PlotUtilities.figure()
        plot_learner_versus_loading_rate_and_number(learner_trials)
        fudge_x_low = 10
        fudge_x_high = 2
        fudge_y = 1.5
        plt.ylim([min_time/fudge_y,max_time*fudge_y])
        plt.xlim([1/fudge_x_low,max(curve_numbers)*fudge_x_high])
        plt.yscale('log')
        plt.xscale('log')        
        PlotUtilities.legend(loc="upper left",frameon=True)
        PlotUtilities.savefig(fig,  base_name + "_all_trials.png")
        # plot the slopes
        fig = PlotUtilities.figure()
        plot_learner_slope_versus_loading_rate(learner_trials)
        PlotUtilities.legend(loc="lower right",frameon=True)
        PlotUtilities.savefig(fig, base_name + "_slopes.png")

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


def autolabel(rects,label_func=lambda i,r: str(r.get_height())):
    """
    Attach a text label above each bar displaying its height

    Args:
        rects: return from ax.bar
        label_func: takes a rect and its index, returs the label
    Returns:
        nothing, but sets text labels
    """
    ax = plt.gca()
    for i,rect in enumerate(rects):
        height = rect.get_height()
        text = label_func(i,rect)
        ax.text(rect.get_x() + rect.get_width()/2., 1.2*height,
                text,ha='center', va='bottom')

def plot_learner_prediction_time_comparison(learners):
    """
    plots the asympotic slope of the learners

    Args:
        learners: list of learners
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
    rects = ax.bar(ind , plot_y, width, color='b',alpha=0.2,linewidth=0,
                   yerr=plot_y_error,log=True,
                   error_kw=dict(ecolor='k',linewidth=2,capsize=10))
    # add some text for labels, title and axes ticks
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(labels)
    xlim = plt.xlim()
    fudge = width/4
    xlim = [min(xlim)-fudge,max(xlim) + fudge]
    plt.xlim(xlim)
    formatted_with_errors = [pretty_exp_with_error(r.get_height(),e) \
                             for r,e in zip(rects,plot_y_error)]
    # add labels with the errors for each
    label_func = lambda i,r : formatted_with_errors[i]
    autolabel(rects,label_func)
    PlotUtilities.lazyLabel("Event finding method",
                            "Points classified per second","")
    
        
def plot_learner_slope_versus_loading_rate(learner_trials):
    """
    Makes a plot of the (slope of runtime versus number of curves) versus
    loading rate

    Args:
        learner_trials: a single learner object
    Returns:
        nothing, makes a pretty plot
    """
    inf = timing_info(learner_trials)
    params,params_std = get_loading_rate_slopes_and_errors(inf)
    # the slope is the time per force extension curve (less an offset; get that
    # per loading rate
    velocities = inf.velocities
    x,xerr = _timing_plot_pts_and_pts_error(inf)
    coeffs,coeffs_err,x_pred_plot,y_pred = get_linear_runtime(x,params)
    # fit a linear model to the runtime
    fudge = 1.75
    plt.errorbar(x=x,xerr=xerr,y=params,yerr=params_std,fmt='ro')
    slope = coeffs[0]
    lower_label = r"{:.2f}$ c_0 N \leq $".format(fudge)
    upper_label = r"$\leq \frac{c_0}{" + "{:.2f}".format(fudge) + "} N$"
    label_timing = lower_label + "T(N)"  + upper_label
    style_timing = dict(color='b',linestyle='--')
    plt.plot(x_pred_plot,y_pred/fudge,label=label_timing,**style_timing)
    plt.plot(x_pred_plot,y_pred*fudge,**style_timing)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    PlotUtilities.lazyLabel("N, Points per curve",
                            "Runtime per curve (s)",
                            "Runtime per curve, T(N), is $\Theta(N)$")



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
    to_ret = r"$" +  sig_fig + r"\cdot 10^{" + sign  + exponent + r"}$" 
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


if __name__ == "__main__":
    run()
