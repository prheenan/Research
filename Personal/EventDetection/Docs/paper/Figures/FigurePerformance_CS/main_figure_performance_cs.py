# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle,os

import svgutils.compose as sc
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import \
    Learning,InputOutput,Plotting,Analysis
from Research.Personal.EventDetection.Util.Learning import  learning_curve

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
    cat= np.concatenate([previous,new])
    to_update_min = np.min(cat)
    if (floor is not None):
        to_update_min = np.max([to_update_min,floor])
    to_update_max = np.max(cat)
    return [to_update_min,to_update_max]

def run(base="./"):
    """
    
    """
    data_file = base + "data/Scores.pkl"
    trials = CheckpointUtilities.lazy_load(data_file)
    out_base = base
    best_fold = []
    lim_load_max = [0,0]
    lim_force_max = [0,0]
    count_max = 0
    for l in trials:
        ruptures_valid_true,ruptures_valid_pred = \
                Learning.get_true_and_predicted_ruptures_per_param(l)
        ret  = [metrics(true,pred) \
                for true,pred in zip(ruptures_valid_true,ruptures_valid_pred)]
        coeffs = [r[0][-1] for r in ret]
        # get the best coefficient
        best_param_idx = np.argmax(coeffs)
        lim_force = ret[best_param_idx][1]
        lim_load = ret[best_param_idx][2]
        counts = ret[best_param_idx][3]
        count_max = max(counts,count_max)
        lim_force_max = update_limits(lim_force,lim_force_max)
        lim_load_max = update_limits(lim_load,lim_load_max,floor=1e-1)
        x_values = l.param_values()
        name = l.description.lower()
        true,pred =  ruptures_valid_true[best_param_idx],\
                     ruptures_valid_pred[best_param_idx]
        best_fold.append([x_values,name,true,pred,best_param_idx])
    # plot the best fold for each
    out_names = []
    for i,(x,name,true,pred,best_param_idx) in enumerate(best_fold):
        # XXX move to checkpointing?
        valid_scores = trials[i]._scores_by_params(train=False)
        distance_scores = Learning.event_distance_f_score(valid_scores)
        best_param_score = distance_scores[best_param_idx]
        # plot everything 
        out_learner_base = "{:s}{:s}".format(out_base,name)
        fig = PlotUtilities.figure(figsize=(8,6))
        Plotting.distance_f_score_plot(best_param_score)
        PlotUtilities.savefig(fig,out_learner_base + "_f.png")
        use_legend = (i == 0)
        fig = PlotUtilities.figure(figsize=(8,6))
        Plotting.rupture_plot(true,pred,use_legend=use_legend,
                              lim_plot_load=lim_load_max,
                              lim_plot_force=lim_force_max,
                              count_limit=[0.5,count_max*2])
        final_out_path = out_learner_base + ".svg"
        PlotUtilities.savefig(fig,final_out_path)
        out_names.append(final_out_path)
    data_panels = [sc.Panel(sc.SVG(f)) for f in out_names]
    sc.Figure("30cm", "37cm", 
              *(data_panels)).\
        tile(1,len(data_panels)).save(out_base + "landscape.svg")
    for f in out_names:
        os.remove(f) 



if __name__ == "__main__":
    run()
