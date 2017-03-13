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


from Research.Personal.EventDetection.Util.Plotting \
    import algorithm_colors,algorithm_markers,algorithm_linestyles


class plotting_metrics:
    def __init__(self,l,ret):
        ruptures_valid_true,ruptures_valid_pred = \
                Learning.get_true_and_predicted_ruptures_per_param(l)
        coeffs = [r[0][-1] for r in ret]
        # get the best coefficient
        self.best_param_idx = np.argmax(coeffs)
        self.lim_force = ret[self.best_param_idx][1]
        self.lim_load = ret[self.best_param_idx][2]
        self.counts = ret[self.best_param_idx][3]
        self.x_values = l.param_values()
        self.name = l.description.lower()
        self.true,self.pred =  ruptures_valid_true[self.best_param_idx],\
                               ruptures_valid_pred[self.best_param_idx]
        self.valid_scores = \
            l._scores_by_params(train=False)[self.best_param_idx]
    def true_and_pred_distances(self,floor_is_max=True):
        kwargs = dict(floor_is_max = floor_is_max)
        to_true = Learning.\
                  event_distance_distribution([self.valid_scores],
                                              to_true=True,**kwargs)[0]
        to_pred = Learning.\
                  event_distance_distribution([self.valid_scores],
                                              to_true=False,**kwargs)[0]
        return to_true,to_pred
    def distance_limit(self,**kwargs):
        true,pred = self.true_and_pred_distances(**kwargs)
        safe_max = lambda x: max(x) if len(x) > 0 else 0
        safe_min = lambda x: min(x) if len(x) > 0 else np.inf
        max_dist = max([safe_max(true),safe_max(pred)])
        min_dist = min([safe_min(true),safe_min(pred)])
        tmp_limits = [min_dist,max_dist]
        return tmp_limits

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


def run(base="./"):
    """
    
    """
    out_base = base
    data_file = base + "data/Scores.pkl"
    force=False
    metric_list = CheckpointUtilities.getCheckpoint(base + "cache.pkl",
                                                    get_metric_list,force,
                                                    data_file)
    lim_load_max = [0,0]
    lim_force_max = [0,0]
    distance_limits = [0,0]
    count_max = 0
    distance_mins = np.inf
    # get the plotting limits
    for m in metric_list:
        lim_force_max = update_limits(m.lim_force,lim_force_max)
        lim_load_max = update_limits(m.lim_load,lim_load_max,floor=1e-1)        
        count_max = max(m.counts,count_max)
        distance_limits = update_limits(m.distance_limit(),distance_limits)
        distance_mins = min(distance_mins,min(m.distance_limit()))
    distance_limits = [distance_mins,max(distance_limits)]
    # POST: have limits...
    # plot the best fold for each
    out_names = []
    colors_pred =  algorithm_colors()
    n_bins = 50
    for i,m in enumerate(metric_list):
        x,name,true,pred = m.x_values,m.name,m.true,m.pred
        best_param_idx = m.best_param_idx
        out_learner_base = "{:s}{:s}".format(out_base,name)
        color_pred =  colors_pred[i]
        # get the distance information we'll need
        to_true,to_pred = m.true_and_pred_distances()
        limit = m.distance_limit()
        log_limit = np.log10(limit)
        bins = np.logspace(*log_limit,num=n_bins)
        fig = PlotUtilities.figure(figsize=(8,8))
        # define the styles for the histograms
        common_style_hist = dict(alpha=0.3,linewidth=0)
        label_true_dist_hist = r"d$_{\mathrm{p}\rightarrow\mathrm{t}}$"
        label_pred_dist_hist = r"d$_{\mathrm{t}\rightarrow\mathrm{p}}$"
        color_true = 'g'
        style_true = dict(color=color_true,label=label_true_dist_hist,
                          **common_style_hist)
        style_pred = dict(color=color_pred,label=label_pred_dist_hist,
                          **common_style_hist)
        distance_histogram = dict(to_true=to_true,to_pred=to_pred,
                                  distance_limits=distance_limits,
                                  bins=bins,style_true=style_true,
                                  style_pred=style_pred)
        Plotting.histogram_event_distribution(**distance_histogram)
        PlotUtilities.savefig(fig,out_learner_base + "dist.png")
        # plot the metric plot
        use_legend = (i == 0)
        fig = PlotUtilities.figure(figsize=(16,8))
        Plotting.rupture_plot(true,pred,use_legend=use_legend,
                              lim_plot_load=lim_load_max,
                              lim_plot_force=lim_force_max,
                              color_pred=color_pred,
                              count_limit=[0.5,count_max*2],
                              distance_histogram=distance_histogram)
        final_out_path = out_learner_base + ".svg"
        PlotUtilities.savefig(fig,final_out_path)
        out_names.append(final_out_path)
    data_panels = [sc.Panel(sc.SVG(f)) for f in out_names]
    sc.Figure("32cm", "41cm", 
              *(data_panels)).\
        tile(1,len(data_panels)).save(out_base + "landscape.svg")
    for f in out_names:
        os.remove(f) 



if __name__ == "__main__":
    run()
