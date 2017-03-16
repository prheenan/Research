# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from Research.Personal.EventDetection.Util.Learning import  learning_curve

from Research.Personal.EventDetection.Util import \
    Learning,InputOutput,Analysis

class coeffs:
    def __init__(self,bc_loading,bc_rupture,bc_2d,median_true,median_pred,
                 q_true,q_pred,name):
        self.name = name
        self.bc_loading = bc_loading
        self.bc_rupture = bc_rupture
        self.bc_2d = bc_2d
        self.median_true=median_true
        self.median_pred=median_pred
        self.q_true=q_true
        self.q_pred=q_pred

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
    def to_true_and_pred_distances(self,floor_is_max=True):
        kwargs = dict(floor_is_max = floor_is_max)
        to_true = Learning.\
                  event_distance_distribution([self.valid_scores],
                                              to_true=True,**kwargs)[0]
        to_pred = Learning.\
                  event_distance_distribution([self.valid_scores],
                                              to_true=False,**kwargs)[0]
        return to_true,to_pred
    def n_true(self):
        f_true = lambda x: x.n_true()
        n_true = Learning.lambda_distribution([self.valid_scores],f_true)[0]
        return n_true
    def n_pred(self):
        f_pred = lambda x: x.n_pred()
        n_pred = Learning.lambda_distribution([self.valid_scores],f_pred)[0]
        return n_pred
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
    def distance_limit(self,**kwargs):
        true,pred = self.to_true_and_pred_distances(**kwargs)
        safe_max = lambda x: max(x[np.where(x>0)]) if len(x) > 0 else -1
        safe_min = lambda x: min(x[np.where(x>0)]) if len(x) > 0 else np.inf
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
        q = 75
        if (len(to_true) > 0):
            median_true = np.median(to_true)
            median_pred  = np.median(to_pred)
            q_true  = np.percentile(to_true,q)
            q_pred = np.percentile(to_pred,q)
        else: 
            median_true,median_pred  = -1,-1
            q_true,q_pred  = -1,-1
        return coeffs(*tmp,median_true=median_true,median_pred=median_pred,
                      q_true=q_pred,q_pred=q_pred,name=self.name)

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

def get_metric_list(data_file):
    trials = CheckpointUtilities.lazy_load(data_file)
    best_fold = []
    metric_list = []
    for l in trials:
        ruptures_valid_true,ruptures_valid_pred = \
                Learning.get_true_and_predicted_ruptures_per_param(l)
        ret  = [metrics(true,pred) \
                for true,pred in zip(ruptures_valid_true,ruptures_valid_pred)]
        tmp = plotting_metrics(l,ret)
        metric_list.append(tmp)
    return metric_list    
