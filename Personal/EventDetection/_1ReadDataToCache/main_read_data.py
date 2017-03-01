# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os

sys.path.append("../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util import Learning
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import Plotting,InputOutput
        


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    cache_directory = "./cache/"
    force = False
    # limit (per category)
    limit = 10
    n_folds = 3
    n_tuning_points = 10
    debug_directory = "./debug_no_event/"
    learners_kwargs = dict(n_points_no_event=n_tuning_points,
                           n_points_fovea=n_tuning_points,
                           n_points_wavelet=n_tuning_points)
    positives_directory = InputOutput.get_positives_directory()
    positive_categories = Learning.get_categories(positives_directory)
    # for each category, predict where events are
    file_name_cache = "{:s}Scores.pkl".format(cache_directory)
    learners = CheckpointUtilities.\
               getCheckpoint(file_name_cache,Learning.get_cached_folds,force,
                             positive_categories,
                             force,cache_directory,limit,n_folds,
                             learners_kwargs=learners_kwargs)
    num_to_plot = 2
    # XXX looking at the worst of the best for the first learner (no event)
    learner = learners[0]
    valid_scores = learner._scores_by_params(train=False)
    x_values = learner.param_values()
    x_tmp,score_tmp,error_tmp = Learning.median_dist_metric(x_values,
                                                            valid_scores)
    # get the lowest median distance ('best case')
    best_x = x_tmp[np.argmin(score_tmp)]
    # find what that means in the real values, if something was invalid
    # (assumes no duplicated params...)
    best_param_idx = np.argmin(np.abs(best_x-x_values))
    # get the corresponding validation folds
    folds = [f for f in learner.validation_folds[best_param_idx]]
    # get all the scores in the folds for the best parameters
    scores = [score for f in folds for score in f.scores]
    print(scores)
    # get all the distances
    distances = [s.minimum_distance_median() for s in scores]
    distances_idx_where_none = np.where(distances is None)[0]
    # get the worst (largest) distances where we arent none
    # XXX note: None is smaller than everything, seems like, so argsort is OK
    sort_idx_high_to_low = np.argsort(distances)[::-1]
    worst_n_idx =  sort_idx_high_to_low[:num_to_plot]
    file_names = [scores[i].source_file + scores[i].name 
                  for i in worst_n_idx]
    # load the worst n back into memory
    # redo the prediction for the worst N, saving to the debug directory
    for l in learners:
        # XXX determine where things went wrong (load/look at specific examples)
        # plot everything
        Plotting.plot_individual_learner(cache_directory,l)



if __name__ == "__main__":
    run()
