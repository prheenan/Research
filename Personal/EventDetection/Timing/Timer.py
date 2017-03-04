# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,random,multiprocessing
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
        self.list_of_time_trials = list_of_time_trials
        self.loading_rates = loading_rates
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
        for i,trial in enumerate(learner_trials.list_of_time_trials):
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
    for c in categories:
        data = Learning.category_read(c,force=True,
                                      cache_directory=cache_directory,
                                      limit=max(curve_numbers))
        c.set_data(data)  
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
    positive_categories = Learning.\
        get_categories(positives_directory=positives_directory)
    curve_numbers = [1,2,5,10,20,35,50,100,150,200]
    cache_dir = "../_1ReadDataToCache/cache/"
    force = False
    times = CheckpointUtilities.getCheckpoint(cache_dir + "all.pkl",
                                              cache_all_learners,force,
                                              learners,positive_categories,
                                              curve_numbers,
                                              cache_dir,force)
    max_time = max([l.max_time_trial() for l in times])
    min_time = min([l.min_time_trial() for l in times])
    for learner_trials in times:
        # plot the timing stuff 
        fig = PlotUtilities.figure()
        plot_learner_versus_loading_rate_and_number(learner_trials)
        fudge = 2
        plt.ylim([min_time/fudge,max_time*fudge])
        plt.xlim([1/fudge,max(curve_numbers)*fudge])
        plt.yscale('log')
        plt.xscale('log')        
        PlotUtilities.legend(loc="lower right",frameon=True)
        PlotUtilities.savefig(fig,learner_trials.learner.description + "_t.png")
        # plot the slopes
        fig = PlotUtilities.figure()
        plot_learner_slope_versus_loading_rate(learner_trials)
        PlotUtilities.legend(loc="lower right",frameon=True)
        PlotUtilities.savefig(fig,learner_trials.learner.description + "_s.png")

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
    coeffs = []
    for num,mean in zip(inf.nums,inf.means):
        coeffs.append(GenUtilities.GenFit(x=num,y=mean))
    # the slope is the time per force extension curve (less an offset; get that
    # per loading rate
    velocities = inf.velocities
    x,xerr = _timing_plot_pts_and_pts_error(inf,round_to_one_decimal=False)
    params = [c[0][0] for c in coeffs]
    params_std = [c[1][0] for c in coeffs]
    plt.errorbar(x=x,xerr=xerr,y=params,yerr=params_std,fmt='ro')
    PlotUtilities.lazyLabel("Thousands of points per curve",
                            "Runtime per curve","")



def _timing_plot_pts_and_pts_error(inf,round_to_one_decimal):
    """
    gets the average number of points per curve in a given loading rate
    and the error (1 standard deviaiton)

    Args:
        inf: the timing_info object
    Returns:
        tuple of <mean number of points per curve, stdev of points per curve>
    """
    n_deci = lambda x: int(np.floor(np.log10(abs(x))))
    rounded = lambda x :int(np.round(x,-n_deci(x)))
    n = inf.pts_per.size
    pts,pts_err = np.zeros(shape=n),np.zeros(shape=n)
    for i,(pts_tmp,xerr_tmp) in enumerate(zip(inf.pts_per,inf.pts_std)):
        pts[i] = pts_tmp/1000
        pts_err[i] = xerr_tmp/1000
    if (not round_to_one_decimal):
        pass
    else:
        # POST: need  to round
        pts = np.array([rounded(p) for p in pts])
        pts_err = np.array([rounded(e) for e in pts_err if e >0])
    return pts,pts_err

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
              dict(color='k',marker='v',linestyle='-.')]
    inf = timing_info(learner_trials)
    pts,xerr = _timing_plot_pts_and_pts_error(inf,True)
    for i,(num,mean,yerr,vel) in enumerate(zip(inf.nums,inf.means,inf.stdevs,
                                               inf.velocities)):
        style = styles[i % len(styles)]
        velocity_label = r"v={:4d}nm/s".format(int(vel))
        number_label = r"N={:d}$\pm${:d}".\
                       format(pts[i],xerr[i])
        label = "{:s}\n({:s})".format(velocity_label,number_label)
        plt.errorbar(x=num,y=mean,yerr=yerr,label=label,**style)
    title = "Runtime verus loading rate and number of curves\n" + \
           "(N, kilopoints/curve, in parenthesis) "
    PlotUtilities.lazyLabel("Number of Force-Extension Curves","Time",title)


if __name__ == "__main__":
    run()
