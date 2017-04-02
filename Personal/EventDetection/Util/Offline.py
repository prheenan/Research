# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from Research.Personal.EventDetection.Util.Learning import  learning_curve
from GeneralUtil.python import CheckpointUtilities,GenUtilities
from Research.Personal.EventDetection.Util import \
    Learning,InputOutput,Analysis

class coeffs:
    def __init__(self,bc_loading,bc_rupture,bc_2d,cat_median,cat_q,
                 cat_relative_median,cat_relative_q,name,q):
        self.name = name
        self.q = q
        self.bc_loading = bc_loading
        self.bc_rupture = bc_rupture
        self.bc_2d = bc_2d
        self.cat_median = cat_median
        self.cat_q = cat_q
        # relative ones!
        self.cat_relative_median=cat_relative_median
        self.cat_relative_q=cat_relative_q

class plotting_metrics:
    def __init__(self,l,ret):
        ruptures_valid_true,ruptures_valid_pred = \
                Learning.get_true_and_predicted_ruptures_per_param(l)
        coeffs = [r[0][-1] for r in ret]
        coeffs_safe_idx = np.where(np.isfinite(coeffs))[0]
        n = len(coeffs)
        idx_arr = np.arange(0,n,step=1)
        coeffs_safe = [coeffs[i] for i in coeffs_safe_idx]
        idx_arr_safe = idx_arr[coeffs_safe_idx]
        best_idx_rel = int(np.round(np.argmax(coeffs_safe)))
        best_idx = idx_arr_safe[best_idx_rel]
        # get the best coefficient
        self._bc_coeffs = coeffs
        self.best_param_idx = best_idx
        self.lim_force = ret[self.best_param_idx][1]
        self.lim_load = ret[self.best_param_idx][2]
        self.counts = ret[self.best_param_idx][3]
        self.x_values = l.param_values()
        self.name = l.description.lower()
        self.true,self.pred =  ruptures_valid_true[self.best_param_idx],\
                               ruptures_valid_pred[self.best_param_idx]
        self.valid_scores = \
            l._scores_by_params(train=False)[self.best_param_idx]
        self.train_scores = \
            l._scores_by_params(train=True)[self.best_param_idx]
    def to_true_and_pred_distances(self,floor_is_max=True):
        kwargs = dict(floor_is_max = floor_is_max)
        to_true = Learning.\
                  event_distance_distribution([self.valid_scores],
                                              to_true=True,**kwargs)[0]
        to_pred = Learning.\
                  event_distance_distribution([self.valid_scores],
                                              to_true=False,**kwargs)[0]
        return to_true,to_pred
    def _lambda(self,f):
        ret= Learning.lambda_distribution([self.valid_scores],f)[0]        
        return ret 
    def max_x_distances_true_pred(self):
        f_max_pred = lambda x: [x.max_displacement() for _ in range(x.n_pred())]
        f_max_true = lambda x: [x.max_displacement() for _ in range(x.n_true())]
        true = np.concatenate(self._lambda(f_max_true))
        pred = np.concatenate(self._lambda(f_max_pred))
        return true,pred
    def n_true(self):
        f_true = lambda x: x.n_true()
        return self._lambda(f_true)
    def n_pred(self):
        f_pred = lambda x: x.n_pred()
        return self._lambda(f_pred)
    def recall(self):
        true = self.n_true()
        pred = self.n_pred()
        true_pos = np.minimum(pred,true)
        return sum(true_pos)/sum(true)
    def precision(self):
        true = self.n_true()
        pred = self.n_pred()
        true_pos = np.minimum(pred,true)
        precision = sum(true_pos)/sum(pred)
        return precision
    def distance_limit(self,relative=False,**kwargs):
        true,pred = self.to_true_and_pred_distances(**kwargs)
        safe_max = lambda x: max(x[np.where(x>0)]) if len(x) > 0 else -1
        safe_min = lambda x: min(x[np.where(x>0)]) if len(x) > 0 else np.inf   
        if (relative):
            max_x_true,max_x_pred =  self.max_x_distances_true_pred()
            true = true/max_x_pred
            pred = pred/max_x_true
        max_dist = max([safe_max(true),safe_max(pred)])
        min_dist = min([safe_min(true),safe_min(pred)])
        tmp_limits = [min_dist,max_dist]

        return tmp_limits
    def coefficients(self):
        ruptures_true,loading_true = \
            Learning.get_rupture_in_pN_and_loading_in_pN_per_s(self.true)
        ruptures_pred,loading_pred = \
            Learning.get_rupture_in_pN_and_loading_in_pN_per_s(self.pred)
        lim_force,bins_rupture,lim_load,bins_load = \
            Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                    loading_true,loading_pred)
        tmp = Analysis.bc_coeffs_load_force_2d(loading_true,loading_pred,
                                               bins_load,ruptures_true,
                                               ruptures_pred,bins_rupture)
        to_true,to_pred = self.to_true_and_pred_distances()
        max_x_true,max_x_pred = self.max_x_distances_true_pred()
        q=85
        cat_median,cat_q,cat_relative_median,cat_relative_q = \
                relative_and_absolute_median_and_q(to_true,to_pred,max_x_true,
                                                   max_x_pred,q=85)
        return coeffs(*tmp,cat_median=cat_median,cat_q=cat_q,name=self.name,q=q,
                      cat_relative_median=cat_relative_median,
                      cat_relative_q=cat_relative_q)

def relative_and_absolute_median_and_q(to_true,to_pred,max_x_true,max_x_pred,
                                       q=85,**kwargs):
    to_true_relative = to_true/max_x_pred
    to_pred_relative = to_pred/max_x_true
    q = 85
    if (len(to_true) > 0):
        cat = np.concatenate([to_true,to_pred])
        cat_rel = np.concatenate([to_true_relative,to_pred_relative])
    else:
        cat = to_true
        cat_rel = to_true_relative
    cat_median = np.median(cat)
    cat_q = np.percentile(cat,q)
    # get the relative metrics
    cat_relative_median = np.median(cat_rel)
    cat_relative_q = np.percentile(cat_rel,85)
    return cat_median,cat_q,cat_relative_median,cat_relative_q

def metrics(true,pred):
    """
    returns useful metrics on a given true and predicted rupture objects

    Args:
        true/pred: the 
    Returns:
        tuple of metric coefficies, limit of force, limits of loading rate,
        and max of (number of true and predicted)
    """
    ruptures_true,loading_true = Learning.\
                                 get_rupture_in_pN_and_loading_in_pN_per_s(true)
    ruptures_pred,loading_pred = Learning.\
                                 get_rupture_in_pN_and_loading_in_pN_per_s(pred)
    lim_force,bins_rupture,lim_load,bins_load = \
        Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                loading_true,loading_pred)
    coeffs = Analysis.\
        bc_coeffs_load_force_2d(loading_true,loading_pred,bins_load,
                                ruptures_true,ruptures_pred,bins_rupture)
    counts = max(len(ruptures_true),len(ruptures_pred))
    return coeffs,lim_force,lim_load,counts

def update_limits(previous,new,floor=None):
    """
    given a previous bounds arrayand a new one, returns an updated to the 
    max of both

    Args:
         previous/new: the two to update from
         floor: if provided, minimum
    Returns:
         updated oubnds
    """
    cat_min= [np.min(previous),np.min(new)]
    cat_max= [np.max(previous),np.max(new)]
    to_update_min = np.min(cat_min)
    if (floor is not None):
        to_update_min = np.max([to_update_min,floor])
    to_update_max = np.max(cat_max)
    return [to_update_min,to_update_max]
def best_metric_from_learner(l):
    """
    returns the best metric from the learner l

    Args:
        l: learning_curve object
    Returns:
        *best* plotting_metrics
    """
    ruptures_valid_true,ruptures_valid_pred = \
        Learning.get_true_and_predicted_ruptures_per_param(l)
    ret  = [metrics(true,pred) \
            for true,pred in zip(ruptures_valid_true,ruptures_valid_pred)]
    return plotting_metrics(l,ret)

def get_metric_list(data_file):
    """
    gets a list of the *best* tuning parameter found from the learners in
    the pickle'd data_file

    Args:
        data_file: the pkl file with all the learnered
    Returns:
        lit of plotting_metric objects (for compression of tuning results)
    """
    trials = CheckpointUtilities.lazy_load(data_file)
    best_fold = []
    metric_list = []
    for l in trials:
        metric_list.append(best_metric_from_learner(l))
    return metric_list    
