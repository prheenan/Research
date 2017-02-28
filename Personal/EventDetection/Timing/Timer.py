# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,random
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
        return [np.mean(x) for x in self.times]
    def std_time_for_fixed_number(self):
        """
        Returns:
            standard deviation of times to complete each of self.num_curves
        """    
        return [np.std(x) for x in self.times]
    def total_number_of_points_per_curves(self):
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

def time_single(func,data):
    """
    time a a single predicton of a set of data, per:
stackoverflow.com/questions/7370801/measure-time-elapsed-in-python/25823885#2582388

    Args:
        func: to call, just takes in the FEC
        data: what to use
    Returns:
        time, in seconds, that it takes to run. 
    """
    start = timer()
    for d in data:
        func(d)
    elapsed_time =  timer() - start
    return elapsed_time

def get_all_times(learner,data,list_of_curve_numbers,trials_per_curve_set=5):
    """
    gets the time_trials for a single loading rate (or whatever)

    Args:
         see cache_all_learners, except a single category's data
    Returns:
        a single time_trials object
    """
    times_all_trials = []
    sizes_all_trials = []
    num_curves_all_trials = []
    for l in list_of_curve_numbers:
        # dont do trials that we cant actually time
        if (l > len(data)):
            continue
        # determine the data set we will use for this one
        data_tmp = data[:l]
        times = []
        sizes = [d.Force.size for d in data_tmp]
        for t in range(trials_per_curve_set):
            times.append(time_single(learner.func_to_call,data_tmp))
        times_all_trials.append(times)
        sizes_all_trials.append(sizes)
        num_curves_all_trials.append(len(data_tmp))
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
    for c in categories:
        data = c.data
        trials.append(get_all_times(learner,data,curve_numbers,**kwargs))
    loading_rates = [c.velocity_nm_s for c in categories]
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
        data = Learning.category_read(c,force=force,
                                      cache_directory=cache_directory,
                                      limit=max(curve_numbers))
        c.set_data(data)  
    # get all the trials for all the learners        
    for l in learners:
        t = CheckpointUtilities.getCheckpoint(l.description,single_learner,
                                              force,l,curve_numbers,categories,
                                              **kwargs)
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
    curve_numbers = [1,2,5,10,20,50]
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
        fig = PlotUtilities.figure()
        plot_single_learner(learner_trials)
        fudge = 2
        plt.ylim([min_time/fudge,max_time*fudge])
        plt.yscale('log')
        PlotUtilities.legend(loc="lower right",frameon=True)
        PlotUtilities.savefig(fig,learner_trials.learner.description + ".png")
        
def plot_single_learner(learner_trials):    
    styles = [dict(color='r',marker='x',linestyle='--'),
              dict(color='b',marker='o',linestyle='-'),
              dict(color='k',marker='v',linestyle='-.')]
    for i,loading_rate_trial in enumerate(learner_trials.list_of_time_trials):
        style = styles[i % len(styles)]
        num_curves = loading_rate_trial.num_curves
        mean = loading_rate_trial.mean_time_for_fixed_number()
        std = loading_rate_trial.std_time_for_fixed_number()
        pts_per_curve = loading_rate_trial.average_number_of_points_per_curve()
        decimal_places = int(np.floor(np.log10(abs(pts_per_curve))))
        round_pts_per_curve = int(np.round(pts_per_curve,-decimal_places))
        velocity = learner_trials.loading_rates[i]
        velocity_label = r"v={:4d}nm/s".format(velocity)
        number_label = r"<Points per curve>={:d}".format(round_pts_per_curve)
        label = "{:s} ({:s})".format(velocity_label,number_label)
        plt.errorbar(x=num_curves,y=mean,yerr=std,label=label,**style)
    PlotUtilities.lazyLabel("Number of Force-Extension Curves","Time","")

if __name__ == "__main__":
    run()
