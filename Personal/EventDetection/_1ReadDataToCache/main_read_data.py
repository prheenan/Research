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
from Research.Personal.EventDetection.Util import Plotting

class ForceExtensionCategory:
    def __init__(self,number,directory,sample,velocity_nm_s,has_events):
        self.category_number = number
        self.directory = directory  
        self.velocity_nm_s = velocity_nm_s
        self.sample = sample
        self.has_events = has_events
        self.data = None
        self.scores = None
    def set_scores(self,scores):
        self.scores = scores
    def set_data(self,data):
        """
        sets the pointer to the list of TimeSepForce objects for this category
        
        Args:
            data: list of TimeSepForce objects
        Returns:
            nothing
        """
        self.data = data 

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    network = FEC_Util.default_data_root()
    base_directory = network + "/4Patrick/CuratedData/Masters_CSCI/"
    positives_directory = base_directory + "Positive/650nm-4x-bio/csv/"
    negatives_directory = "XXX TODO"
    cache_directory = "./cache/"
    # tuple of <relative directory,sample,velocity> for FEC with events
    positive_meta = \
      [[positives_directory + "1000-nanometers-per-second/","650nm DNA",1000],
       [positives_directory + "500-nanometers-per-second/","650nm DNA",500], 
       [positives_directory + "100-nanometers-per-second/","650nm DNA",100]]
    # create objects to represent our data categories
    positive_categories = [ForceExtensionCategory(i,*r,has_events=True) 
                           for i,r in enumerate(positive_meta)]
    force = True
    debug_plots = True
    # limit (per category)
    limit = 30
    n_folds = 3
    # for each category, predict where events are
    file_name_cache = "{:s}Scores.pkl".format(cache_directory)
    learners = CheckpointUtilities.\
               getCheckpoint(file_name_cache,Learning.get_cached_folds,force,
                             positive_categories,
                             force,cache_directory,limit,n_folds)
    for l in learners:
        Plotting.plot_individual_learner(cache_directory,l)



if __name__ == "__main__":
    run()
