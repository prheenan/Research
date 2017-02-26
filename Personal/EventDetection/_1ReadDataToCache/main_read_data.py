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
    def __init__(self,directory,velocity_nm_s,sample,has_events):
        self.directory = directory  
        self.velocity_nm_s = velocity_nm_s
        self.sample = sample
        self.has_events = has_events
        self.data = None
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
      [[positives_directory + "1000-nanometers-per-second/","650nm DNA",100],
       [positives_directory + "500-nanometers-per-second/","650nm DNA",500], 
       [positives_directory + "100-nanometers-per-second/","650nm DNA",1000]]
    # create objects to represent our data categories
    positive_categories = [ForceExtensionCategory(*r,has_events=True) 
                           for r in positive_meta]
    force = True
    # limit (per category)
    limit = 2
    # get the positive events
    try:
        InputOutput.set_and_cache_category_data(positive_categories,
                                                cache_directory=cache_directory,
                                                force=force,limit=limit)
    except OSError:
        # just read in the files that live here XXX just for debugging
        all_files = GenUtilities.getAllFiles(cache_directory,ext=".pkl")
        single = positive_categories[0]
        all_files = [CheckpointUtilities.getCheckpoint(f,None,False)
                     for f in all_files]
        single.set_data(all_files)
    thresh = 1e-2                                
    splits, prediction_info,scores = [],[],[]
    # for each category, predict where events are
    for i,category in enumerate(positive_categories):
        if (category.data is None):
            continue  
        for j,example in enumerate(category.data):
            example_split = Analysis.\
                zero_and_split_force_extension_curve(example)
            m_func =Detector.adhesion_function_for_split_fec(example_split)
            info = Detector._predict_helper(example_split,threshold=thresh,
                                            condition_function=m_func)
            score = Scoring.get_scoring_info(example_split,info.event_idx)
            prediction_info.append(info)
            splits.append(example_split)
            scores.append(score)
    rupture_true = [r for s in scores for r in s.ruptures_true ]
    rupture_predicted = [r for s in scores for r in s.ruptures_predicted ]
    fig = PlotUtilities.figure(figsize=(8,8))
    Plotting.plot_predicted_and_true_ruptures(rupture_true,rupture_predicted)
    PlotUtilities.savefig(fig,cache_directory + "rupture_loading.png")
    for i,(example_split,info,score) in \
        enumerate(zip(splits,prediction_info,scores)):
        true,pred = score.ruptures_true,score.ruptures_predicted
        out_file_path = cache_directory + "{:d}".format(i)
        fig = PlotUtilities.figure(figsize=(8,20))
        Plotting.plot_autocorrelation(example_split)
        PlotUtilities.savefig(fig,out_file_path + "auto.png")   
        # XXX fix threshhold
        fig = PlotUtilities.figure(figsize=(8,12))    
        Plotting.plot_prediction_info(example_split,info)
        PlotUtilities.savefig(fig,out_file_path + "info.png")
    # get the negative events
    # XXX 


if __name__ == "__main__":
    run()
