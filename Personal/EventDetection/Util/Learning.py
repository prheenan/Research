# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,traceback

from Research.Personal.EventDetection.Util import Analysis,InputOutput,Scoring
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from sklearn.model_selection import StratifiedKFold
import multiprocessing


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
    def param_values(self):
        # XXX assume 1-D search, only one parameters per list
        return np.array([l.values()[0] for l in self.list_of_params])
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

def rupture_objects(scores,get_true):
    """
    get the rupture objects associated with the scores

    Args: 
         scores: see event_distance_distribution
         get_true: if true, gets the *true* rupture objects associated with the
    Returns: 
         array of rupture objects, one per parameter
    """
    if (get_true):
        func_tmp = lambda x: [v.ruptures_true for v in x]
    else: 
        func_tmp = lambda x: [v.ruptures_predicted for v in x]
    # need to concatenate everything
    func_fold = lambda *args,**kwargs: np.concatenate(func_tmp(*args,**kwargs))
    return _walk_scores(scores,func_fold=func_fold,
                        func_param=np.concatenate,func_top=np.array)



def limits_and_bins_force_and_load(ruptures_true,ruptures_pred,
                                   loading_true,loading_pred,n=20):
    """
    Return a 4-tuple of limit,bins  for rupture force and loading rate

    Args:
        <x>_<true/pred> : llist of true/predicted x
        n: number of bins
    Returns:
       limits force,bins force,limits loaidng,bins loading
    """
    double_f = lambda f,*args: f([f(x) for x in args if len(x) > 0])
    # determine the limits on the rupture force
    min_y = double_f(min,ruptures_pred,ruptures_true)
    max_y = double_f(max,ruptures_pred,ruptures_true)
    lim_force = [min_y,max_y]
    # determine the limits on the loading rate
    safe = lambda x: [x[i] for i in np.where(np.array(x)>0)[0]]
    min_x = double_f(min,safe(loading_pred),safe(loading_true))
    max_x = double_f(max,safe(loading_pred),safe(loading_true))
    lim_load = [min_x,max_x]
    bins_rupture= np.linspace(*lim_force,num=n)
    min_y = max(min(lim_load),1e-2)
    logy = np.log10([min_y,max(lim_load)])
    bins_load = np.logspace(*logy,num=n)
    return lim_force,bins_rupture,lim_load,bins_load
    

def get_rupture_in_pN_and_loading_in_pN_per_s(objs):
    """
    Args:
        objs: see _plot_rupture_objecs
    Returns:
        tuple of <rupture force in pN, loading rate in pN>
    """
    to_pN = lambda x: x * 1e12
    rupture_forces_pN = np.array([to_pN(obj.rupture_force) for obj in objs])
    loading_rate_pN_per_s = np.array([to_pN(obj.loading_rate) for obj in objs])
    return rupture_forces_pN,loading_rate_pN_per_s

def get_true_and_predicted_ruptures_per_param(learner):
    """
    gets the truee and preicted rupture objects for the *validation* folds
    of each learner object 

    Args:
         learner: the learner_curve obect to use
    Returns:
         tuple of validation true, predicted ruptures
    """
    train_scores = learner._scores_by_params(train=True)
    valid_scores = learner._scores_by_params(train=False)
    # get the validation ruptures (both truee and predicted)
    ruptures_valid_true = rupture_objects(valid_scores,get_true=True)
    ruptures_valid_pred = rupture_objects(valid_scores,get_true=False)
    return ruptures_valid_true,ruptures_valid_pred


def event_distance_distribution(scores,**kwargs):
    """
    gets the distribution of distances at each paramter value

    Args:
         scores: learner._scores_by_params object
         **kwargs: passed to minimum_distance_distribution
    Returns:
         concatenates distributions at each parameter value
    """
    func_fold = lambda x: \
        np.concatenate([v.minimum_distance_distribution(**kwargs) for v in x])
    return _walk_scores(scores,func_fold = func_fold,
                        func_param=np.concatenate,func_top=np.array)

def f_score_dist(v):
    """
    returns the distance f score for the given score object v

    Args:
        v: to scoore
    Returns;
        f score, 0 to 1, higher is better. 
    """
    kw = dict(floor_is_max=True)
    dist_to_true = v.minimum_distance_distribution(to_true=True,**kw)
    dist_to_pred = v.minimum_distance_distribution(to_true=False,**kw)
    max_x = v.max_x
    # get the averages ? XXX
    if (len(dist_to_true) != 0):
        average_to_true = np.median(dist_to_true)
    else:
        average_to_true = 0
    if (len(dist_to_pred) != 0):
        average_to_pred = np.median(dist_to_pred)
    else:
        average_to_pred = 0
    # defining precision and recall in a distance-analogy sense
    precision_dist = average_to_true
    recall_dist = average_to_pred
    f_score = \
      1-(2/max_x) * (precision_dist * recall_dist)/(precision_dist+recall_dist)
    return f_score

def event_distance_f_score(scores):
    """
    returns the *distance* f score for each curve

    Args:
        scores: see fold_number_events_off
    """
    func_fold = lambda x: [f_score_dist(v) for v in x]
    return _walk_scores(scores,func_fold = func_fold,
                        func_param=np.concatenate,func_top=np.array)
    

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
    relative_missed = np.sum(missed_list)/np.sum(true_list)
    to_ret = relative_missed
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

def median_dist_metric(x_values,scores,**kwargs):
    """
    function for safely getting the median metric

    Args:
        x_values: the parameters
        scores: see safe_scores
        **kwargs: passed to minimum_distance_median...
    Returns:
        see valid_scores_erors_and_params
    """
    score_func_pred = lambda x: median_dist_per_param(x,**kwargs)
    error_func_pred = lambda x: stdev_dist_per_param(x,**kwargs)
    kw_pred = dict(score_func=score_func_pred,error_func=error_func_pred)
    x,dist,dist_std = valid_scores_erors_and_params(x_values,scores,**kw_pred)
    return x,dist,dist_std

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
        
def single_example_info_and_score(func,example,**kwargs):
    """
    Returns a list of learning_curve objects

    Args:
        func: a method taking in example and **kwargs and returning event
        indices
        
        example: TimeSepForce object
        **kwargs: for func
    Returns:
        tuple of <scoring object, folding_meta object>
    """
    meta = example.Meta
    # get the predicted event index
    event_idx = func(example,**kwargs)
    example_split = Analysis.zero_and_split_force_extension_curve(example)
    # get the score
    score = Scoring.get_scoring_info(example_split,event_idx)
    return score,fold_meta(meta)

def run_functor(functor):
    """
    Given a no-argument functor, run it and return its result. We can 
    use this with multiprocessing.map and map it over a list of job 
    functors to do them.

    Handles getting more than multiprocessing's pitiful exception output
    
    Args:
        functor: a no-argument function (probably a lambda)
    Returns:
        whatever functor does, possibly raising an exception
    """
    try:
        # This is where you do your actual work
        return functor()
    except:
        # Put all exception text into an exception and raise that
        err_string = "".join(traceback.format_exception(*sys.exc_info()))
        raise Exception(err_string)    
    
def single_example_multiproc(args):
    """
    multiprocesing interface to single_example_info_and_score

    Args:
        tuple of arguments to single_example_info_and_score
    Returns:
         see single_example_info_and_score
    """
    func,example,dict_kwargs = args
    return single_example_info_and_score(func,example,**dict_kwargs)        
        
class multiprocessing_functor(object):
    def __init__(self):
        pass
    def __call__(self,*args):
        return single_example_multiproc(*args)
            
def single_fold_score(fold_data,func,kwargs,pool):
    """
    Gets the fold object for a single set of data (ie: a single fold)
    for a fixed parameter set

    Args:
        fold_data: set of TimeSepForce objects to use
        func: to call, should return a list of predicted event starts
        kwargs: dict, fixed parameters to pass to func
        pool: if not none, a Multiprocessing pool for parallel processing...
    Returns:
        fold object
    """
    scores = []
    info = []
    if (pool is None):
        scores_info = [single_example_info_and_score(func,ex,**kwargs) 
                       for ex in fold_data]
    else:
        # we make a list of functor objects (functions + arguments)
        # that we can then safely run
        functors_args = [ (func,ex,kwargs) for ex in fold_data]
        scores_info = pool.map(multiprocessing_functor(),functors_args)     
    # POST: got the scores an info, somehow...
    scores = [s[0] for s in scores_info]
    info = [s[1] for s in scores_info]
    return fold(kwargs,scores,info)

def folds_for_one_param(data,param,fold_idx,func_to_call,pool):
    folds = []
    folds_valid = []
    for train_idx,test_idx in fold_idx:
        # get the training scores
        fold_data = [data[f] for f in train_idx]
        folds_tmp = single_fold_score(fold_data,func_to_call,param,
                                      pool=pool)
        folds.append(folds_tmp)
        # get the validation scores
        valid_data = [data[f] for f in test_idx]
        valid_fold = single_fold_score(valid_data,func_to_call,param,
                                       pool=pool)
        folds_valid.append(valid_fold)    
    return folds,folds_valid
    
def get_all_folds_for_one_learner(cache_directory,force,learner,data,fold_idx,
                                  pool):
    """
    Gets all the folds for a single learner

    Args:
        cache_directory: base where we cache things
        force: if true, force re-reading of all folds
        learner: learning_curve to use
        data: list of TimeSepForce objects to use
        fold_idx: list of <train,test>; one tuple per fold
        pool: a multiprocessing pool to draw from (if not None)
    Returns:
        tuple of :
         (0) list, one element per paramter. each element is a list of folds
         (1) list, one element per parameter. each element is a list of 
         validation folds
    """
    func_to_call = learner.func_to_call
    params_then_folds,param_validation_fold = [],[]
    for i,param in enumerate(learner.list_of_params):   
        param_val = param.values()[0]
        cache_name = "{:s}_{:s}_param_{:d}_{:.3g}.pkl".\
           format(cache_directory,learner.description,i,param_val)
        ret = CheckpointUtilities.getCheckpoint(cache_name,folds_for_one_param,
                                                force,data,param,fold_idx,
                                                func_to_call,pool=pool)
        folds,folds_valid = ret
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
    

def get_single_learner_folds(cache_directory,force,l,data,fold_idx,pool_size):
    """
    return the training and testing folds for a given learner

    Args:
        cache_directory: where to cache individual parameter folds
        force: if caching should be forced
        l : learner to use
        d : data to use 
        fold_idx: which indices to use for the folds
        pool_size: number of processing to use. If <=1, just uses 1 
        (no parallelism)
    Returns:
        tuple of <training,validation folds>
    """
    if (pool_size <= 1):
        # dont use a multiprocessing pool
        pool = None
    else:
        pool = multiprocessing.Pool(pool_size)
    print("Using {:d} processes for {:s}".format(pool_size,l.description))
    ret = get_all_folds_for_one_learner(cache_directory,force,
                                        l,data,fold_idx,pool=pool)
    list_of_folds,validation_folds = ret
    return  list_of_folds,validation_folds   

def get_cached_folds(categories,force_read,force_learn,
                     cache_directory,limit,n_folds,seed=42,
                     learners=None,pool_size=1):
    """
    caches all the results for every learner after reading in all the data

    Args:
    
        categories: list of velocity-separated data
        force_read: if the csv fiels should be re-read
        force_learn: if the learner objects should be re-read
        cache_directoy: where to put the pkl files
        limit: how many examples to read in 
        n_folds: now many folds to use
        seed: for PRNG
    Returns:
        list, one element per paramter. each element is a list of folds
    """
    # read and update all the categories
    categories = InputOutput.\
                 read_categories(categories,force_read,cache_directory,limit)
    labels_data = [ [i,d] for i,cat in enumerate(categories) for d in cat.data]
    labels = [l[0] for l in labels_data]
    data = [l[1] for l in labels_data]
    # determine the folds to use
    fold_obj = StratifiedKFold(n_splits=n_folds,shuffle=True,random_state=seed)
    # .split returns a generator by default; convert to a list to avoid
    # making it only used for the first fold 
    fold_idx = list(fold_obj.split(X=np.zeros(len(labels)),y=labels))
    if (learners is None):
        learners = get_learners()
    # POST: all data read in. get all the scores for all the learners.
    for l in learners:
        cache_file = cache_directory + "folds_" + l.description + ".pkl"
        tmp = CheckpointUtilities.getCheckpoint(cache_file,
                                                get_single_learner_folds,
                                                force_learn,
                                                cache_directory,force_learn,
                                                l,data=data,fold_idx=fold_idx,
                                                pool_size=pool_size)
        list_of_folds,validation_folds = tmp
        l.set_list_of_folds(list_of_folds)
        l.set_validation_folds(validation_folds)
    return learners
