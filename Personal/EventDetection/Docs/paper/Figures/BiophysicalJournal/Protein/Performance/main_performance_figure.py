# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting,Learning


def coeff_str(i,j,c,coeffs_array,func):
    to_ret = "{:.3g}".format(c)
    if (func(coeffs_array[:,j]) == i):
        to_ret = r"**" + to_ret + "**"
    return to_ret

def write_coeffs_file(out_file,coeffs):
    opt_low = lambda x: np.argmin(x)
    opt_high = lambda x: np.argmax(x)
    funcs_names_values = []
    for c in coeffs:
        q = c.q
        err_str =  (r"Relative event error $P_{" + "{:d}".format(q) + \
                    r"}$ ($\downarrow$)")
        tmp =[ [opt_low,r"Rupture BCC ($\downarrow$)",
                1-c.bc_2d],
               [opt_low,err_str,c.cat_relative_q]]
        funcs_names_values.append(tmp)
    # only get the funcs nad names from the first (redudant to avoid typos
    funcs = [coeff_tmp[0] for coeff_tmp in funcs_names_values[0] ]
    coeff_names = [coeff_tmp[1] for coeff_tmp in funcs_names_values[0] ]
    plot_dict = Plotting.algorithm_title_dict()
    method_names = [plot_dict[c.name] for c in coeffs]
    # get the list of coefficients
    coeffs = [[f[2] for f in coeff_tmp] for coeff_tmp in funcs_names_values ]
    coeffs_array = np.array(coeffs)
    join_str = " | "
    str_v = join_str + join_str.join(coeff_names) 
    for i,(name,c) in enumerate(zip(method_names,coeffs)): 
        str_v += "\n"
        str_v += "{:s}{:s}".format(name,join_str)
        str_v += join_str.join([coeff_str(i,j,c_tmp,coeffs_array,funcs[j]) 
                                for j,c_tmp in enumerate(c)])
    # add the preamble...
    final = "{:s}\n".format(str_v) + r"[{#tbl:Performance}]" + "\n"
    with open(out_file,'w') as f:
        f.write(final)

def tick_style_log(**kwargs):
    PlotUtilities.tom_ticks(plt.gca(),num_major=6,**kwargs)


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    data_file = "../_Data/Scores.pkl"
    metrics = CheckpointUtilities.getCheckpoint("./cache.pkl",
                                                Offline.get_best_metrics,False,
                                                data_file)
    coeffs_compare = [m.coefficients() for m in metrics]
    write_coeffs_file("./coeffs.txt",coeffs_compare)
    fig = PlotUtilities.figure(figsize=(7,4))
    colors = Plotting.algorithm_colors()
    n_rows = 3
    n_cols = 3
    xlim_dist = [1e-5,2]
    xlim_load = [1,1e5]
    xlim_rupture = [-5,300]
    legend_locs = ['upper right','upper left','upper left']
    titles = ["FEATHER","Fovea","Wavelet"]
    for i,m in enumerate(metrics):
        offset = n_rows * i
        # the first column gets the algorithm label; the first row gets the
        # metric label
        kw_tmp = dict(title_kwargs=dict(fontweight='bold',color='b',fontsize=9),
                      legend_kwargs=dict(fontsize=8,handlelength=0.75,
                                         handletextpad=0.25))
        if offset == 0:
            title_dist = "Location error"
            title_load = r"Loading rate (NuG2 + $\mathrm{\alpha}_3$D)"
            title_rupture_force = r"Rupture force (NuG2 + $\mathrm{\alpha}_3$D)"
        else:
            title_dist,title_load,title_rupture_force = "","",""
        # only have an x label on the last row
        last_row = (offset/n_rows == n_rows-1)
        if (last_row):
            xlabel_dist = "Relative Error ($\mathbf{x_k}$)"
            xlabel_load = "Loading Rate (pN/s)"
            xlabel_rupture_force = "F$_R$ (pN)"
        else:
            xlabel_dist, xlabel_load,xlabel_rupture_force = "","",""
        ylabel_dist = \
            (r"$N_{\mathrm{" + "{:s}".format(titles[i]) + "}}$")
        color_pred=colors[i]
        color_true = 'g'
        # get the formatting dictionaries for the various plots 
        distance_histogram_kw= \
            Plotting.event_error_kwargs(m,color_pred=color_pred,
                                        label_bool=False)
        true_style_histogram = Plotting.\
            _histogram_true_style(color_true=color_true,label="true")
        pred_style_histogram = Plotting.\
            _histogram_predicted_style(color_pred=color_pred,label="predicted")
        # get the binning for the rupture force and loading rates
        true,pred = m.true,m.pred        
        ruptures_true,loading_true = \
            Learning.get_rupture_in_pN_and_loading_in_pN_per_s(true)
        ruptures_pred,loading_pred = \
            Learning.get_rupture_in_pN_and_loading_in_pN_per_s(pred)
        _lim_force_plot,_bins_rupture_plot,_lim_load_plot,_bins_load_plot = \
            Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                    loading_true,loading_pred,
                                                    limit=0.02)            
        # # make the 'just the distance' figures
        ax_dist = plt.subplot(n_rows,n_cols,(offset+1))
        Plotting.histogram_event_distribution(use_q_number=True,
                                              **distance_histogram_kw)
        PlotUtilities.lazyLabel(xlabel_dist,ylabel_dist,title_dist,
                                loc=legend_locs[i],legendBgColor='w',
                                frameon=False,**kw_tmp)
        PlotUtilities.ylabel(ylabel_dist,fontweight='bold')
        plt.xlim(xlim_dist)                                   
        if not last_row:
            PlotUtilities.no_x_label(ax_dist)
        tick_style_log()
        # # make the loading rate histogram      
        ax_load = plt.subplot(n_rows,n_cols,(offset+2))
        Plotting.loading_rate_histogram(pred,bins=_bins_load_plot,
                                        **pred_style_histogram)
        Plotting.loading_rate_histogram(true,bins=_bins_load_plot,
                                        **true_style_histogram)
        plt.xscale('log')
        plt.yscale('log')
        PlotUtilities.lazyLabel(xlabel_load,"",title_load,
                                legendBgColor='w',
                                loc='upper left',frameon=False,**kw_tmp)
        plt.xlim(xlim_load)               
        if not last_row:
            PlotUtilities.no_x_label(ax_load)
        PlotUtilities.no_y_label(ax_load)
        tick_style_log()
        # # make the rupture force histogram
        ax_rupture = plt.subplot(n_rows,n_cols,(offset+3))
        Plotting.rupture_force_histogram(pred,bins=_bins_rupture_plot,
                                        **pred_style_histogram)
        Plotting.rupture_force_histogram(true,bins=_bins_rupture_plot,
                                        **true_style_histogram)
        plt.yscale('log')
        PlotUtilities.lazyLabel(xlabel_rupture_force,"",title_rupture_force,
                                useLegend=False,**kw_tmp)       
        plt.xlim(xlim_rupture)                                                 
        if not last_row:
            PlotUtilities.no_x_label(ax_rupture)     
        PlotUtilities.no_y_label(ax_rupture)                
        tick_style_log(change_x=False)           
        # set all the y limits for this row
        axes_counts = [ax_rupture,ax_load,ax_dist]
        max_limit = np.max([r.get_ylim() for r in axes_counts])
        ylim_new = [0.5,max_limit*1.6]
        for r in axes_counts:
            r.set_ylim(ylim_new)
    axis_func = lambda axes: [ax for i,ax in enumerate(axes) if i < 3]
    loc_last_two = [-0.05,1.1]
    locs = [ [-0.25,1.1], loc_last_two,loc_last_two]
    PlotUtilities.label_tom(fig,axis_func=axis_func,loc=locs)
    PlotUtilities.savefig(fig,"./performance.png",
                          subplots_adjust=dict(hspace=0.1,wspace=0.1))
    

if __name__ == "__main__":
    run()
