# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os

sys.path.append("../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util import \
    Analysis,Plotting,InputOutput,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities

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

def debugging_plots(id_string,example_split,info):
    """
    Plots the autocorrelation and prediction information

    Args:
        id_string: 
        example_split: force extension curve to debug
        info: return from predict_helper
    Returns:
        nothing, splits plots out like id_string
    """
    out_file_path =  id_string
    fig = PlotUtilities.figure(figsize=(8,20))
    Plotting.plot_autocorrelation(example_split)
    PlotUtilities.savefig(fig,out_file_path + "auto.png")   
    # XXX fix threshhold
    fig = PlotUtilities.figure(figsize=(8,12))    
    Plotting.plot_prediction_info(example_split,info)
    PlotUtilities.savefig(fig,out_file_path + "info.png")


def category_read(category,force,cache_directory,limit):
    """
    Reads in all the data associated with a category

    Args:
        category: ForceExtensionCategory object
        force: if true, force re-reading
        cache_directory: if force is not true, where to re-read from
        limit: maximum number to re-read
    Returns:
        list of TimeSepForce objects
    """
    try:
        return InputOutput.get_category_data(category,force,cache_directory,
                                             limit)
    except OSError:
        if (category.category_number != 0):
            return []
        # just read in the files that live here XXX just for debugging
        file_names = GenUtilities.getAllFiles(cache_directory,ext="csv.pkl")
        all_files = [CheckpointUtilities.getCheckpoint(f,None,False)
                     for f in file_names]
        return all_files

        
def category_scores(category,force,cache_directory,limit,debug_plots,
                    thresh):
    """
    Gets all the scores associated with a given category

    Args:
        see category_read, except
        thresh: threshold we want
        debug_plots: if true, save plots out
    Returns:
        list scoring objects
    """
    id_data = "{:s}{:.1f}t={:.5f}".format(category.sample,
                                          category.velocity_nm_s,
                                          thresh)
    cache_data = cache_directory + id_data + ".pkl"
    data = CheckpointUtilities.getCheckpoint(cache_data,category_read,
                                             force,category,force,
                                             cache_directory,limit)
    category.data = data
    scores = []
    if (len(category.data)  == 0):
        return scores
    for j,example in enumerate(category.data):
        example_split = Analysis.zero_and_split_force_extension_curve(example)
        m_func =Detector.adhesion_function_for_split_fec(example_split)
        info = Detector._predict_helper(example_split,threshold=thresh,
                                        condition_function=m_func)
        score = Scoring.get_scoring_info(example_split,info.event_idx)
        scores.append(score)
        wave_name = example_split.retract.Meta.Name
        id_string = cache_directory + "db_" + id_data + "_" + \
                    str(j) + wave_name 
        if (debug_plots):
            debugging_plots(id_string,example_split,info)
    return scores

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
      [[positives_directory + "1000-nanometers-per-second/","650nm DNA",100],
       [positives_directory + "500-nanometers-per-second/","650nm DNA",500], 
       [positives_directory + "100-nanometers-per-second/","650nm DNA",1000]]
    # create objects to represent our data categories
    positive_categories = [ForceExtensionCategory(i,*r,has_events=True) 
                           for i,r in enumerate(positive_meta)]
    force = True
    debug_plots = True
    # limit (per category)
    limit = 2
    thresh = 1e-2                                
    # for each category, predict where events are
    for i,category in enumerate(positive_categories):
        kw_scores = (category,force,cache_directory,
                     dict(limit=limit,thresh=thresh,debug_plots=debug_plots))
        file_name_cache = "{:s}{:s}{:.1f}Scores.pkl".\
            format(cache_directory,category.sample,category.velocity_nm_s)
        scores = CheckpointUtilities.\
            getCheckpoint(file_name_cache,category_scores,force,
                          category,force,cache_directory,
                          limit=limit,thresh=thresh,debug_plots=debug_plots)
        category.set_scores(scores)
    scores_flat = [score for cat in positive_categories for score in cat.scores]
    rupture_true = [r for s in scores_flat for r in s.ruptures_true ]
    rupture_predicted = [r for s in scores_flat for r in s.ruptures_predicted ]
    fig = PlotUtilities.figure(figsize=(8,8))
    Plotting.plot_predicted_and_true_ruptures(rupture_true,rupture_predicted)
    PlotUtilities.savefig(fig,cache_directory + "rupture_loading.png")


if __name__ == "__main__":
    run()
