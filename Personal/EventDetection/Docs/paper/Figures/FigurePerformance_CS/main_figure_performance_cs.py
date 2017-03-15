# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle,os,string

import svgutils.compose as sc
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import \
    Learning,InputOutput,Plotting,Analysis
from Research.Personal.EventDetection.Util.Learning import  learning_curve
import matplotlib.gridspec as gridspec


from Research.Personal.EventDetection.Util.Plotting \
    import algorithm_colors,algorithm_markers,algorithm_linestyles


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

def write_coeffs_file(out_file,coeffs):
    opt_low = lambda x: np.argsort(x)
    opt_high = lambda x: opt_low(x)[::-1]
    funcs_names_values = []
    for c in coeffs:
        tmp =[ [opt_high,"bc2",c.bc_2d],
               [opt_low,"medtrue",c.median_true],
               [opt_low,"medpred",c.median_pred],
               [opt_low,"qtrue",c.q_true],
               [opt_low,"qpred",c.q_pred]]
        funcs_names_values.append(tmp)
    # only get the funcs nad names from the first (redudant to avoid typos
    funcs = [coeff_tmp[0] for coeff_tmp in funcs_names_values[0] ]
    coeff_names = [coeff_tmp[1] for coeff_tmp in funcs_names_values[0] ]
    method_names = [c.name for c in coeffs]
    # get the list of coefficients
    coeffs = [[f[2] for f in coeff_tmp] for coeff_tmp in funcs_names_values ]
    join_str = " & "
    str_v = join_str + join_str.join(coeff_names) + "\e\\hline \n"
    for name,c in zip(method_names,coeffs): 
        str_v += "{:s}{:s}".format(name,join_str)
        str_v += join_str.join(["{:.3g}".format(c_tmp) for c_tmp in c])
        str_v += "\\e\n"
    with open(out_file,'w') as f:
        f.write(str_v)

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
    coeffs_compare = []
    # get the plotting limits
    for m in metric_list:
        lim_force_max = update_limits(m.lim_force,lim_force_max)
        lim_load_max = update_limits(m.lim_load,lim_load_max,floor=1e-1)        
        count_max = max(m.counts,count_max)
        distance_limits = update_limits(m.distance_limit(),distance_limits)
        distance_mins = min(distance_mins,min(m.distance_limit()))
        coeffs_compare.append(m.coefficients())
    write_coeffs_file(out_base + "coeffs.txt",coeffs_compare)
    distance_limits = [distance_mins,max(distance_limits)]
    # POST: have limits...
    # plot the best fold for each
    out_names = []
    colors_pred =  algorithm_colors()
    n_bins = 50
    for m in metric_list:
        precision = m.precision()
        recall = m.recall()
        data = [precision,recall]
        bins = np.linspace(0,1)
        colors = ['g','r']
        labels = ["Precision","Recall"]
        ind = np.arange(len(data)) 
        style_precision = dict(color='g',alpha=0.3,hatch="//",label="Precision")
        style_recall = dict(color='r',alpha=0.3,label="Recall")
        width = 0.5
        rects = plt.bar(ind,data,color=colors,width=width)
        ax = plt.gca()
        ax.set_xticks(ind + width / 2)
        ax.set_xticklabels(labels)
        y_func=lambda i,r: "{:.3f}".format(r.get_height()/2)
        PlotUtilities.autolabel(rects,y_func=y_func)
        PlotUtilities.lazyLabel("Metric Value",
                                "Number of Force-Extension Curves","",
                                frameon=True)
        plt.xlim([-0.1,1+2*width])
        plt.ylim([0,1])
        #XXX todo
    # make a giant figure, 3 rows (one per algorithm)
    fig = PlotUtilities.figure(figsize=(16,22))
    entire_figure = gridspec.GridSpec(3,1)
    for i,m in enumerate(metric_list):
        x,name,true,pred = m.x_values,m.name,m.true,m.pred
        best_param_idx = m.best_param_idx
        out_learner_base = "{:s}{:s}".format(out_base,name)
        color_pred =  colors_pred[i]
        # get the distance information we'll need
        to_true,to_pred = m.to_true_and_pred_distances()
        limit = m.distance_limit()
        log_limit = np.log10(limit)
        bins = np.logspace(*log_limit,num=n_bins)
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
        gs = gridspec.GridSpecFromSubplotSpec(2, 3, width_ratios=[2,2,1],
                                              height_ratios=[2,2,1],
                                              subplot_spec=entire_figure[i],
                                              wspace=0.25,hspace=0.2)
        # plot the metric plot
        use_legend = (i == 0)
        Plotting.rupture_plot(true,pred,use_legend=use_legend,
                              lim_plot_load=lim_load_max,
                              lim_plot_force=lim_force_max,
                              color_pred=color_pred,
                              count_limit=[0.5,count_max*2],
                              distance_histogram=distance_histogram,gs=gs,
                              fig=fig)
    # individual plot labels
    n_subplots = 5
    n_categories = len(metric_list)
    letters =  string.lowercase[:n_categories]
    letters = [ ["({:s}{:d})".format(s,n+1) for n in range(n_subplots)]
                 for s in letters]
    flat_letters = [v for list_of_v in letters for v in list_of_v]
    PlotUtilities.label_tom(fig,flat_letters,loc=(-1.1,1.1),fontsize=15)
    final_out_path = out_base + "landscape.pdf"
    PlotUtilities.savefig(fig,final_out_path)




if __name__ == "__main__":
    run()
