# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,random,multiprocessing,re
from timeit import default_timer as timer

sys.path.append("../../../../")
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import Learning,InputOutput,Learners
from Research.Personal.EventDetection.Timing import TimePlot

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
    def __init__(self,learner,list_of_time_trials,loading_rates,curve_numbers):
        self.learner = learner
        # sort the time trials, sorting increasing by number of points
        idx_sort = np.argsort([t.average_number_of_points_per_curve
                               for t in list_of_time_trials])
        self.list_of_time_trials = [list_of_time_trials[i] for i in idx_sort]
        self.loading_rates = [loading_rates[i] for i in idx_sort]
        self.curve_numbers = curve_numbers
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
    return time_trials_by_loading_rate(learner,trials,loading_rates,
                                       curve_numbers)

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
    learners = Learners.get_learners()
    positives_directory = InputOutput.get_positives_directory()
    positive_categories = InputOutput.\
        get_categories(positives_directory=positives_directory,
                       use_simulated=True)
    curve_numbers = [1,2,5,10,30,50,100,150,200]
    cache_dir = "../_1ReadDataToCache/cache/"
    GenUtilities.ensureDirExists(cache_dir)
    force = False
    times = CheckpointUtilities.getCheckpoint(cache_dir + "all_timer.pkl",
                                              cache_all_learners,force,
                                              learners,positive_categories,
                                              curve_numbers,
                                              cache_dir,force)
    out_base = "./out/"
    GenUtilities.ensureDirExists(out_base)
    # sort the times by their loading rates
    max_time = max([l.max_time_trial() for l in times])
    min_time = min([l.min_time_trial() for l in times])
    # plot the Theta(n) coefficient for each
    fig = PlotUtilities.figure()
    TimePlot.plot_learner_prediction_time_comparison(times)
    PlotUtilities.legend(loc="lower right",frameon=True)
    PlotUtilities.savefig(fig,out_base + "compare.png")
    for learner_trials in times:
        base_name = out_base + learner_trials.learner.description
        # plot the timing veruses loading rate and number of points 
        fig = PlotUtilities.figure()
        TimePlot.plot_learner_versus_loading_rate_and_number(learner_trials)
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
        TimePlot.plot_learner_slope_versus_loading_rate(learner_trials)
        PlotUtilities.legend(loc="lower right",frameon=True)
        PlotUtilities.savefig(fig, base_name + "_slopes.png")

if __name__ == "__main__":
    run()
