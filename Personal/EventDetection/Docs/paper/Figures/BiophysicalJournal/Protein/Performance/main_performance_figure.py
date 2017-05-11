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
    fig = PlotUtilities.figure(figsize=(18,12))
    colors = Plotting.algorithm_colors()
    n_rows = 3
    n_cols = 3
    xlim_dist = [1e-5,2]
    xlim_load = [1,1e5]
    xlim_rupture = [-5,300]
    legend_locs = ['upper right','upper left','upper right']
    for i,m in enumerate(metrics):
        offset = n_rows * i
        # the first column gets the algorithm label; the first row gets the
        # metric label
        if offset == 0:
            title_dist = "Location Error"
            title_load = "Loading Rate Histogram"
            title_rupture_force = "Rupture Force Histogram"
        else:
            title_dist,title_load,title_rupture_force = "","",""
        # only have an x label on the last row
        last_row = (offset/n_rows == n_rows-1)
        if (last_row):
            xlabel_dist = "Relative Error (x$_k$)"
            xlabel_load = "Loading Rate (pN/s)"
            xlabel_rupture_force = "Rupture Force (pN)"
        else:
            xlabel_dist, xlabel_load,xlabel_rupture_force = "","",""
        ylabel_dist = "{:s}\nCount".format(m.name.title())
        color_pred=colors[i]
        color_true = 'g'
        # get the formatting dictionaries for the various plots 
        distance_histogram_kw= \
            Plotting.event_error_kwargs(m,color_pred=color_pred)
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
                                                    limit=0.01)            
        # make the 'just the distance' figures
        ax_dist = plt.subplot(n_rows,n_cols,(offset+1))
        Plotting.histogram_event_distribution(use_q_number=True,
                                              **distance_histogram_kw)
        PlotUtilities.lazyLabel(xlabel_dist,ylabel_dist,title_dist,
                                loc=legend_locs[i],legendBgColor='w',
                                frameon=True)      
        plt.xlim(xlim_dist)                                   
        if not last_row:
            PlotUtilities.no_x_label()
        # make the loading rate histogram      
        ax_load = plt.subplot(n_rows,n_cols,(offset+2))
        Plotting.loading_rate_histogram(pred,bins=_bins_load_plot,
                                        **pred_style_histogram)
        Plotting.loading_rate_histogram(true,bins=_bins_load_plot,
                                        **true_style_histogram)
        plt.xscale('log')
        plt.yscale('log')
        PlotUtilities.lazyLabel(xlabel_load,"",title_load,
                                legendBgColor='w',
                                loc='upper left',frameon=True)
        plt.xlim(xlim_load)               
        if not last_row:
            PlotUtilities.no_x_label()
        PlotUtilities.no_y_label()
        # make the rupture force histogram
        ax_rupture = plt.subplot(n_rows,n_cols,(offset+3))
        Plotting.rupture_force_histogram(pred,bins=_bins_rupture_plot,
                                        **pred_style_histogram)
        Plotting.rupture_force_histogram(true,bins=_bins_rupture_plot,
                                        **true_style_histogram)
        plt.yscale('log')
        PlotUtilities.lazyLabel(xlabel_rupture_force,"",title_rupture_force,
                                useLegend=False)       
        plt.xlim(xlim_rupture)                                                 
        if not last_row:
            PlotUtilities.no_x_label()     
        PlotUtilities.no_y_label()                           
        # set all the y limits for this row
        axes_counts = [ax_rupture,ax_load,ax_dist]
        max_limit = np.max([r.get_ylim() for r in axes_counts])
        ylim_new = [0.5,max_limit*5]
        for r in axes_counts:
            r.set_ylim(ylim_new)
    PlotUtilities.savefig(fig,"./performance.png")
    

if __name__ == "__main__":
    run()
