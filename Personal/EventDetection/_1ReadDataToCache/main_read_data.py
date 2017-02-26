# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os
from sklearn.cross_validation import StratifiedKFold

sys.path.append("../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util import \
    Analysis,Plotting,InputOutput,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities

class fold:
    def __init__(self,params,scores,info):
        self.params = params
        self.scores = scores
        self.info = info

class learning_curve:
    def __init__(self,list_of_folds,list_of_params,func_to_call):
        self.list_of_folds = list_of_folds
        self.list_of_params = list_of_params
        self.func_to_call = func_to_call
    def _concatenate_all_scores(self):
        return [score for list_by_params in self.list_of_folds
                for s in list_by_params for score in s.scores]
            

class learners:
    def __init__(self,no_event,fovea,wavelet,fold_info):
        self.no_event = no_event
        self.fovea = fovea
        self.wavelet = wavelet

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

def single_fold_score(category,fold_data,func_to_call,*args,**kwargs):
    """
    XXX add in actual function support, remove extra ars
    """
    cache_directory = "./cache/"
    scores = []
    params = []
    info = []
    for j,example in enumerate(fold_data):
        sample = category[j].sample,
        velocity = category[j].velocity_nm_s
        id_data = "{:s}{:.1f}p={:s}".format(sample,velocity,str(kwargs))
        example_split = Analysis.zero_and_split_force_extension_curve(example)
        m_func =Detector.adhesion_function_for_split_fec(example_split)
        final_dict = dict(condition_function=m_func,**kwargs)
        pred_info = Detector._predict_helper(example_split,**final_dict)
        score = Scoring.get_scoring_info(example_split,pred_info.event_idx)
        scores.append(score)
        wave_name = example_split.retract.Meta.Name
        id_string = cache_directory + "db_" + id_data + "_" + \
                    str(j) + wave_name 
        debugging_plots(id_string,example_split,pred_info)
        info.append([sample,velocity])
        params.append(kwargs)
    return fold(params,scores,info)


def get_fold_scores(categories,n_folds=2,
                    params_no_event=[dict(threshold=0.01)]):
    labels_data = [ [i,d] for i,cat in enumerate(categories) for d in cat.data]
    labels = [l[0] for l in labels_data]
    data = [l[1] for l in labels_data]
    # determine the folds to use
    fold_idx = StratifiedKFold(labels,n_folds)
    params_then_folds = []
    for p in params_no_event:
        folds = []
        for train_idx,test_idx in fold_idx:
            category_tmp = [categories[labels[f]] for f in train_idx]
            fold_data = [data[f] for f in train_idx]
            folds_tmp = single_fold_score(category_tmp,fold_data,
                                          func_to_call=None,**p)
            folds.append(folds_tmp)
        params_then_folds.append(folds)
    # XXX fix parameter stuff
    return learning_curve(params_then_folds,
                          list_of_params=params_no_event,func_to_call=None)

def get_cached_folds(categories,force,cache_directory,limit,thresh):
    for c in categories:
        data_tmp = category_read(c,force,cache_directory,limit)
        c.set_data(data_tmp)
    # POST: all data read in. get all the scores for all the learners.
    return get_fold_scores(categories)

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
    force = False
    debug_plots = True
    # limit (per category)
    limit = 2
    thresh = 1e-2                                
    # for each category, predict where events are
    file_name_cache = "{:s}Scores.pkl".format(cache_directory)
    learners = CheckpointUtilities.\
               getCheckpoint(file_name_cache,get_cached_folds,force,
                             positive_categories,
                             force,cache_directory,limit,thresh)
    scores = learners._concatenate_all_scores()
    print(len(scores))
    rupture_true = [r for s in scores for r in s.ruptures_true ]
    rupture_predicted = [r for s in scores for r in s.ruptures_predicted ]
    fig = PlotUtilities.figure(figsize=(8,8))
    Plotting.plot_predicted_and_true_ruptures(rupture_true,rupture_predicted)
    PlotUtilities.savefig(fig,cache_directory + "rupture_loading.png")


if __name__ == "__main__":
    run()
