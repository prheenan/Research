# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from Research.Personal.EventDetection.Util import Analysis,InputOutput,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from Research.Personal.EventDetection.OtherMethods.Roduit2012_OpenFovea import \
    fovea
from Research.Personal.EventDetection.OtherMethods.Numpy_Wavelets import \
    wavelet_predictor
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
    def _scores_by_params(self,train=True):
        fold_list = self.list_of_folds if train else self.validation_folds
        scores_by_params = [ [[s for s in fold.scores]
                              for fold in folds_by_param]
                             for folds_by_param in fold_list]
        return scores_by_params

def _walk_scores(scores,
                 func_score=lambda x : x,
                 func_fold =lambda x : x,
                 func_param=lambda x : x,
                 func_top  = lambda x:x):
    """
    function for 'easily' walking through list of scores

    Args:
        func_<x>: applied at the level of x (e.g. func_score is for a single
        score, func_fold is a list of scores meaning a fold, etc)
    Returns:
         result of the chained functions
    """
    return  func_top([ func_param([ func_fold([func_score(s) for s in scores])
                          for scores in by_param])
                       for by_param in scores])


def safe_scores(scores,value_func=lambda x: x,eval_func=lambda x:x):
    """
    function for getting possibly None values and evaluating them

    Args:
        scores: list of scores
        value:func: takes a score, gives a value
        eval_func: for evaluating the non-None values
    Returns:
         result of the chained functions
    """
    raw = [value_func(s) for s in scores]
    safe = [r for r in raw if r is not None]
    if (len(safe) > 0):
        return eval_func(safe)
    else:
        return None

def safe_median(scores):
    """
    function for safely evaluating the median

    Args:
        scores: see safe_scores
    Returns:
         result of the chained functions
    """
    return safe_scores(scores,eval_func=np.median)

def median_dist_per_param(scores,**kwargs):
    """
    function for safely getting the median of the scores we want

    Args:
        scores: see safe_scores
        **kwargs: passed to minimum_distance_median
    Returns:
        median of the minimum distance to an event, per paramter across folds
        (1-D arrray)
    """
    score_func = lambda x: x.minimum_distance_median(**kwargs)
    func_fold = lambda x: safe_scores(x,value_func=score_func,
                                      eval_func=np.median)
    return _walk_scores(scores,func_fold =func_fold,
                        func_param=safe_median,func_top=np.array)

def stdev_dist_per_param(scores,**kwargs):
    """
    function for safely getting the median of the scores we want

    Args:
        scores: see safe_scores
        kwargs: see median_dist_per_param
    Returns:
        stdev of the minimum distance to an event, per paramter across folds
        (1-D arrray)
    """
    score_func = lambda x: x.minimum_distance_distribution(**kwargs)
    eval_func = lambda x: np.std(np.concatenate(x))
    func_fold = lambda x: safe_scores(x,value_func=score_func,
                                      eval_func=eval_func)
    return _walk_scores(scores,func_fold =func_fold,
                        func_param=safe_median,func_top=np.array)

def fold_number_events_off(scores):
    """
    see number_events_off_per_param

    Args:
         scores: see number_events_off_per_param
         learning_curve._scores_by_params
    returns:
         sum of missed events divided by sum of true events
    """
    true_pred = [x.n_true_and_predicted_events() for x in scores]
    true_list = [t[0] for t in true_pred]
    pred_list = [t[1] for t in true_pred]
    missed_list = [abs(true-pred) for true,pred in zip(true_list,pred_list)]
    relative_missed = np.array(missed_list)/np.array(true_list)
    to_ret = np.median(relative_missed)
    return to_ret

def number_events_off_per_param(params,scores):
    """
    gets the (relative) number of events we were off by:
    (1) gettting the predicted and true number of events in a fold
    (2) getting rel=missed - true
    (3) taking the mean and stdev of rel across all folds

    Args:
         params: the x value to use
         scores: the scorer object to use, formatted like 
         learning_curve._scores_by_params
    returns;
         tuple of <valid params, valid scores, valie errors>
    """
    cat_median = lambda x: safe_median(np.concatenate(x))
    cat_std = lambda x: np.std(np.concatenate(x))
    score_func = lambda s : _walk_scores(s,func_fold=fold_number_events_off,
                                         func_param=safe_median,
                                         func_top=np.array)
    error_func = lambda s:  _walk_scores(s,func_fold=fold_number_events_off,
                                         func_param=np.std,
                                         func_top=np.array)
    kw = dict(score_func=score_func,error_func=error_func)
    return valid_scores_erors_and_params(params,scores,**kw)

def valid_scores_erors_and_params(params,scores,score_func,error_func):
    """
    given a function for getting scores and errors, finds where the results
    are valid, selecting the x,y, and erro there

    Args:
        params: what the x values are 
        scores: the scores we will search using score/error_func
        <score/error>_func: functions giving the scores and errors of scores
        at each value of params. If undefined, then is None
    Returns:
        tuple of <valid x, valid score, valid score error>
    """
    dist = score_func(scores)
    dist_std = error_func(scores)
    valid_func = lambda x: (~np.equal(x,None))
    good_idx_func = lambda train,valid : np.where( valid_func(train) & \
                                                   valid_func(valid))
    good_idx = good_idx_func(dist,dist_std)
    return params[good_idx],dist[good_idx],dist_std[good_idx]


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
    scores = []
    info = []
    for j,example in enumerate(fold_data):
        meta = example.Meta
        # get the predicted event index
        event_idx = func(example,**kwargs)
        example_split = Analysis.zero_and_split_force_extension_curve(example)
        # get the score
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
        # done with all the folds for this parameter; save them out
        params_then_folds.append(folds)
        param_validation_fold.append(folds_valid)
    return params_then_folds,param_validation_fold

def _get_single_curve(name,tuple_v,func):
    """
    Returns a single learning curve object

    Args:
        name: the name of the curvess
        tuple_v: tuple like <function to call, list of parameters>
            
        func: takes in list of single parameters, returns a list of kwargs  
        dicts for func
    Returns:
        learning_curve object
    """
    return learning_curve(name,tuple_v[0],func(tuple_v[1]))
    
def get_learners(n_points_no_event=5,n_points_fovea=5,n_points_wavelet=5):
    """
    Returns a list of learning_curve objects

    Args:
        n_points_no_event: number of points for varying the no event portion of 
        things
            
        n_points_fovea: number of points to use on fovea
    Returns:
        list of learning curves
    """
    # make the no event example
    no_event_func = lambda arg_list: [dict(threshold=t) for t in arg_list]
    no_event_tuple = [Detector.predict,np.logspace(-3.5,-1.5,endpoint=True,
                                                   base=10,
                                                   num=n_points_no_event)]
    no_event_curve = _get_single_curve("No Event",no_event_tuple,no_event_func)                                
    # make the fovea example
    fovea_func = lambda arg_list: [dict(weight=w) for w in arg_list]
    fovea_tuple = [fovea.predict,np.linspace(0.01,0.5,endpoint=True,
                                             num=n_points_fovea)]
    fovea_curve = _get_single_curve("Open Fovea",fovea_tuple,fovea_func)                                   
    # make the CWT example
    cwt_func = lambda arg_list: [dict(min_snr=w) for w in arg_list]
    cwt_tuple = [wavelet_predictor.predict,
                 np.linspace(start=20,stop=150,num=n_points_wavelet)]
    wavelet_curve = _get_single_curve("Wavelet transform",cwt_tuple,cwt_func)   
    return [no_event_curve,fovea_curve,wavelet_curve]

def get_single_learner_folds(l,data,fold_idx):
    """
    return the training and testing folds for a given learner

    Args:
        l : learner to use
        d : data to use 
        fold_idx: which indices to use for the folds
    Returns:
        tuple of <training,validation folds>
    """
    list_of_folds,validation_folds = get_all_folds_for_one_learner(l,data,
                                                                   fold_idx)
    return  list_of_folds,validation_folds                                                                                 

    
def get_cached_folds(categories,force,cache_directory,limit,n_folds,seed=42,
                     learners_kwargs=dict()):
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
    learners = get_learners(**learners_kwargs)
    # POST: all data read in. get all the scores for all the learners.
    for l in learners:
        cache_file = cache_directory + "folds_" + l.description + ".pkl"
        tmp = CheckpointUtilities.getCheckpoint(cache_file,
                                                get_single_learner_folds,force,
                                                l,data,fold_idx)
        list_of_folds,validation_folds = tmp
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
    return positive_categories
