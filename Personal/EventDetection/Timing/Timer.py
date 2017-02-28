# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,random

from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import Learning

class time_trials:
    def __init__(self,times,num_curves,sizes):
        self.times = times
        self.num_curves = num_curves
        self.fec_num_points = fec_num_points

class time_trials_by_loading_rate:
    def __init__(self,learner,list_of_time_trials,loading_rates):
        self.learner = learner
        self.list_of_time_trials = list_of_time_trials
        self.loading_rates = loading_rates


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
    for l in list_of_curve_numbers:
        shuffled_data = random.shuffle(data)
        if (l < len(shuffled_data)):
            continue
        data_tmp = shuffled_data[:l]
        times = []
        sizes = [d.Force.size for d in data_tmp]
        for t in trials_per_curve_set:
            times.append(time_single(learner.func_to_call,data_tmp))
        times_all_trials.append(times)
        sizes_all_trials.append(sizes)
        num_curves_all_trials.append(len(data_tmp))
    return time_trial(times_all_trials,num_curves_all_trials,sizes_all_trials)

def single_learner(learner,curve_numbers,categories,**kwargs):
    """
    gets the time_trials_by_loading_rate for a single object

    Args:
         see cache_all_learners, except a single learner
    Returns:
        a single time_trials_by_loading_rate
    """
    for c in categories:
        data = c.data
        trials.append(get_all_times(learner,data,curve_numbers,**kwargs))
    loading_rates = [c.velocity_nm_s for c in categories]
    return time_trials_by_loading_rate(learner,trials,loading_rates)

def cache_all_learners(learners,categories,curves_numbers,**kwargs):
    """
    caches and returns the timing results for all the learners

    Args:
        learners: list ofLearning.learnin_curve object to use
        categories: list of Learninrg.ForceExtensionCtegories to use
        curve_numbers: how many curves to use (from each category;
        they are all used separately)
        **kwargs: passed to curve_numbers
    Returns:
        list of time_trials_by_loading_rate objects
    """
    times = []
    for l in learners:
        t = CheckpointUtilities.getCheckpoint(l.description,
                                              l,curve_numbers,categories,
                                              **kwargs)
        times.append(t)
    return t
        
        
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    learners = Learning.get_learners()
    positive_categories = Learning.get_categories(positives_directory)
    curve_numbers = [1,2,5,10,20,50,100,200]
    for c in positive_categories:
        Learning.category_read(c,force=False,cache_directory="./cache/",
                               limit=max(curve_numbers))
    times = cache_all_learners(learners,positive_categories,curves_numbers)
    print(times[0].list_of_time_trials[0].times)

if __name__ == "__main__":
    run()
