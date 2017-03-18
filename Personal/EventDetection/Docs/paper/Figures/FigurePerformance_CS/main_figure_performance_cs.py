# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle,os,string

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting,Analysis
from Research.Personal.EventDetection.Util.Offline import plotting_metrics
import matplotlib.gridspec as gridspec

from Research.Personal.EventDetection.Util.Plotting \
    import algorithm_colors,algorithm_markers,algorithm_linestyles

def coeff_str(i,j,c,coeffs_array,func):
    to_ret = "{:.3g}".format(c)
    if (func(coeffs_array[:,j]) == i):
        to_ret = r"\textbf{" + to_ret + "}"
    return to_ret

def write_coeffs_file(out_file,coeffs):
    opt_low = lambda x: np.argmin(x)
    opt_high = lambda x: np.argmax(x)
    funcs_names_values = []
    for c in coeffs:
        tmp =[ [opt_high,r"bc2 ($\uparrow$)",c.bc_2d],
               [opt_low,r"median [nm] ($\downarrow$)",c.cat_median*1e9],
               [opt_low,r"q [nm]($\downarrow$)",c.cat_q*1e9],
               [opt_low,r"relative median ($\downarrow$)",
                c.cat_relative_median],
               [opt_low,r"relative q ($\downarrow$)",c.cat_relative_q]]
        funcs_names_values.append(tmp)
    # only get the funcs nad names from the first (redudant to avoid typos
    funcs = [coeff_tmp[0] for coeff_tmp in funcs_names_values[0] ]
    coeff_names = [coeff_tmp[1] for coeff_tmp in funcs_names_values[0] ]
    method_names = [c.name for c in coeffs]
    # get the list of coefficients
    coeffs = [[f[2] for f in coeff_tmp] for coeff_tmp in funcs_names_values ]
    coeffs_array = np.array(coeffs)
    join_str = " & "
    str_v = join_str + join_str.join(coeff_names) + "\e\\hline \n"
    for i,(name,c) in enumerate(zip(method_names,coeffs)): 
        str_v += "{:s}{:s}".format(name,join_str)
        str_v += join_str.join([coeff_str(i,j,c_tmp,coeffs_array,funcs[j]) 
                                for j,c_tmp in enumerate(c)])
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
                                                    Offline.get_metric_list,
                                                    force,data_file)
    lim_load_max = [0,0]
    lim_force_max = [0,0]
    distance_limits = [0,0]
    count_max = 0
    distance_mins = np.inf
    coeffs_compare = []
    # get the plotting limits
    update_limits = Offline.update_limits
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
    fig = PlotUtilities.figure(figsize=(16,18))
    entire_figure = gridspec.GridSpec(3,1)
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 4
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
        label_pred_dist_hist = r"d$_{\mathrm{p}\rightarrow\mathrm{t}}$"
        label_true_dist_hist = r"d$_{\mathrm{t}\rightarrow\mathrm{p}}$"
        color_true = 'g'
        style_true = dict(color=color_true,label=label_true_dist_hist,
                          **common_style_hist)
        style_pred = dict(color=color_pred,label=label_pred_dist_hist,
                          **common_style_hist)
        use_legend = (i == 0)
        xlabel_histogram = "Distance [m]" if (i == len(metric_list)-1) else ""
        distance_histogram = dict(to_true=to_true,to_pred=to_pred,
                                  distance_limits=distance_limits,
                                  bins=bins,style_true=style_true,
                                  style_pred=style_pred,
                                  xlabel=xlabel_histogram)
        gs = gridspec.GridSpecFromSubplotSpec(2, 3, width_ratios=[2,2,1],
                                              height_ratios=[2,1],
                                              subplot_spec=entire_figure[i],
                                              wspace=0.4,hspace=0.5)
        # plot the metric plot
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
    PlotUtilities.label_tom(fig,flat_letters,loc=(-0.22,1.1),fontsize=18)
    final_out_path = out_base + "landscape.pdf"
    PlotUtilities.savefig(fig,final_out_path,
                          subplots_adjust=dict(left=0.10,
                                               hspace=0.2,wspace=0.2,top=0.95))




if __name__ == "__main__":
    run()
