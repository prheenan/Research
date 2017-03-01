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
    positives_directory = InputOutput.get_positives_directory()
    cache_directory = "./cache/"
    positive_categories = Learning.get_categories(positives_directory)
    force = False
    debug_plots = True
    # limit (per category)
    limit = 5
    n_folds = 3
    n_tuning_points = 10
    learners_kwargs = dict(n_points_no_event=n_tuning_points,
                           n_points_fovea=n_tuning_points,
                           n_points_wavelet=n_tuning_points)
    # for each category, predict where events are
    file_name_cache = "{:s}Scores.pkl".format(cache_directory)
    learners = CheckpointUtilities.\
               getCheckpoint(file_name_cache,Learning.get_cached_folds,force,
                             positive_categories,
                             force,cache_directory,limit,n_folds,
                             learners_kwargs=learners_kwargs)
    for l in learners:
        Plotting.plot_individual_learner(cache_directory,l)



if __name__ == "__main__":
    run()
