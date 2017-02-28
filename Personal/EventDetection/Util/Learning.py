# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from Research.Personal.EventDetection.Util import Analysis,InputOutput,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from sklearn.cross_validation import StratifiedKFold

class fold_meta:
    def __init__(self,meta):
        self.velocity = meta.Velocity
        self.name = meta.Name
        self.source_file = meta.SourceFile

class fold:
    def __init__(self,param,scores,info):
        """
        stores a single 'fold' (ie: fixed parameter value, scores and meta
        for each FEC)

        Args:
            param: parameter used
            scores: list of score object, one per FEC in thefold
            info: meta information about the FEC objects
        """
        self.param = param
        self.scores = scores
        self.info = info

class learning_curve:
    def __init__(self,description,func_to_call,list_of_params):
        """
        stores the folds and parameters associated with a function and list of 
        learning parameters

        Args:
            description: short name
            func_to_call: what function to call for this data
            list of paramters: list of dictionaries (Each of which is 
            passed to func_to_call for the appropriate folds)
        """
        self.description = description
        self.list_of_folds = None
        self.validation_folds = None
        self.list_of_params = list_of_params
        self.func_to_call = func_to_call
    def set_list_of_folds(self,folds):
        self.list_of_folds = folds
    def set_validation_folds(self,folds):
        self.validation_folds = folds
    def _concatenate_all_scores(self):
        fold_score_list = [f for folds in self.list_of_folds for f in folds]
        all_scores = [s for fold in fold_score_list for s in fold.scores]
        return all_scores

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

def single_fold_score(fold_data,func,kwargs):
    """
    Gets the fold object for a single set of data (ie: a single fold)
    for a fixed parameter set

    Args:
        fold_data: set of TimeSepForce objects to use
        func: to call, should return a list of predicted event starts
        **kwargs: fixed parameters to pass to func
    Returns:
        fold object
    """
    cache_directory = "./cache/"
    scores = []
    info = []
    for j,example in enumerate(fold_data):
        meta = example.Meta
        event_idx = func(example,**kwargs)
        example_split = Analysis.zero_and_split_force_extension_curve(example)
        score = Scoring.get_scoring_info(example_split,event_idx)
        info.append(fold_meta(meta))
        scores.append(score)
    return fold(kwargs,scores,info)


def get_all_folds_for_one_learner(learner,data,fold_idx):
    """
    Gets all the folds for a single learner

    Args:
        learner: learning_curve to use
        data: list of TimeSepForce objects to use
        fold_idx: list of <train,test>; one tuple per fold
    Returns:
        tuple of :
         (0) list, one element per paramter. each element is a list of folds
         (1) list, one element per parameter. each element is a list of 
         validation folds
    """
    func_to_call = learner.func_to_call
    params_then_folds,param_validation_fold = [],[]
    for param in learner.list_of_params:
        folds = []
        folds_valid = []
        for train_idx,test_idx in fold_idx:
            # get the training scores
            fold_data = [data[f] for f in train_idx]
            folds_tmp = single_fold_score(fold_data,func_to_call,param)
            folds.append(folds_tmp)
            # get the validation scores
            valid_data = [data[f] for f in test_idx]
            valid_fold = single_fold_score(valid_data,func_to_call,param)
            folds_valid.append(valid_fold)
        params_then_folds.append(folds)
        param_validation_fold.append(folds_valid)
    return params_then_folds,param_validation_fold

def get_learners(n_points_no_event=7):
    """
    Returns a list of learning_curve objects

    Args:
        n_points_no_event: number of points for varying the no event portion of 
        things
    Returns:
        list of lerning curvess
    """
    no_event_args_to_dict = lambda arg_list: \
            [dict(threshold=t) for t in arg_list]
    no_event_tuple = [Detector.predict,np.linspace(1e-4,1e-1,endpoint=True,
                                                   num=n_points_no_event)]
    no_event_curve = learning_curve("No event",no_event_tuple[0],
                                    no_event_args_to_dict(no_event_tuple[1]))
    return [no_event_curve]

def get_cached_folds(categories,force,cache_directory,limit,n_folds,seed=42):
    """
    caches all the results for every learner after reading in all the data

    Args:
        categories: list of velocity-separated data
        force: if the csv fiels should be re-read
        cache_directoy: where to put the pkl files
        limit: how many examples to read in 
        n_folds: now many folds to use
        seed: for PRNG
    Returns:
        list, one element per paramter. each element is a list of folds
    """
    for c in categories:
        data_tmp = category_read(c,force,cache_directory,limit)
        c.set_data(data_tmp)
    labels_data = [ [i,d] for i,cat in enumerate(categories) for d in cat.data]
    labels = [l[0] for l in labels_data]
    data = [l[1] for l in labels_data]
    # determine the folds to use
    fold_idx = StratifiedKFold(labels,n_folds=n_folds,shuffle=True,
                               random_state=seed)
    learners = get_learners()
    # POST: all data read in. get all the scores for all the learners.
    for l in learners:
        list_of_folds,validation_folds = \
            get_all_folds_for_one_learner(l,data,fold_idx)
        l.set_list_of_folds(list_of_folds)
        l.set_validation_folds(validation_folds)
    return learners

def get_categories(positives_directory):
    """
    get all the categories associated with the loading rates we will use

    Args:
        positives_directory: base directory where things live
    Returns:
        list of ForceExtensionCategory
    """
    # tuple of <relative directory,sample,velocity> for FEC with events
    positive_meta = \
      [[positives_directory + "1000-nanometers-per-second/","650nm DNA",1000],
       [positives_directory + "500-nanometers-per-second/","650nm DNA",500], 
       [positives_directory + "100-nanometers-per-second/","650nm DNA",100]]
    # create objects to represent our data categories
    positive_categories = [ForceExtensionCategory(i,*r,has_events=True) 
                           for i,r in enumerate(positive_meta)]
    
