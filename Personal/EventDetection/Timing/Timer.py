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
    def mean_and_stdev_time_for_fixed_number_of_curves(self):
        walk = lambda f : [f(x) for x in self.times]
        return walk(np.mean),walk(np.std)
    def total_number_of_points_per_curves(self):
        return [sum(n) for n in self.fec_num_points]
    

class time_trials_by_loading_rate:
    def __init__(self,learner,list_of_time_trials,loading_rates):
        self.learner = learner
        self.list_of_time_trials = list_of_time_trials
        self.loading_rates = loading_rates

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
    force = True
    times = cache_all_learners(learners,positive_categories,curve_numbers,
                               force=force,cache_directory=cache_dir)
    for learner in times:
        analyze_single_learner(learner)
        
def analyze_single_learner(learner_trials):    
    fig = PlotUtilities.figure()
    for i,loading_rate_trial in enumerate(learner_trials.list_of_time_trials):
        num_curves = loading_rate_trial.num_curves
        mean,std =  \
            loading_rate_trial.mean_and_stdev_time_for_fixed_number_of_curves()
        print(i,mean,std)
        plt.errorbar(x=num_curves,y=mean,yerr=std,fmt='ro-')
    PlotUtilities.lazyLabel("Number of Force-Extension Curves","Time","")
    PlotUtilities.savefig(fig,learner_trials.learner.description + ".png")

if __name__ == "__main__":
    run()
