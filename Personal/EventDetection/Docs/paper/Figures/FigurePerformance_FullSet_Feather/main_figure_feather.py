# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle,os,string

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting

def get_best_metrics(data_file):
    learners = CheckpointUtilities.lazy_load(data_file)
    metrics = [Offline.best_metric_from_learner(l) for l in learners]
    return metrics

def run(in_base="./"):
    """
    
    """
    out_base = "./"
    data_file = in_base + "data/Scores.pkl"
    force=False
    cache_file = out_base + "cache.pkl"
    metrics = CheckpointUtilities.getCheckpoint(cache_file,get_best_metrics,
                                                force,data_file)
    loc_left = (-0.15,1.1)
    loc_top = (-0.15,1.05)
    loc_lower = (-0.15,0.95)
    locs = [loc_top,loc_left,loc_left,loc_lower,loc_lower]
    titles = Plotting.algorithm_title_dict()
    colors = Plotting.algorithm_colors()
    for i,m in enumerate(metrics):
        name = titles[m.name.lower()]
        safe_name = name.replace(" ","")
        color_pred =  colors[i]
        distance_histogram= Plotting.event_error_kwargs(m,color_pred=color_pred)
        true,pred = m.true,m.pred
        # make the 'just the distance' figures
        fig = PlotUtilities.figure((10,6))
        Plotting.histogram_event_distribution(**distance_histogram)
        final_out_dist = "{:s}{:s}_dist_full.pdf".format(out_base,safe_name)
        PlotUtilities.savefig(fig,final_out_dist)
        # make the rupture spectrum figure
        fig = PlotUtilities.figure((12,7))
        final_out_rupture = "{:s}{:s}_rupture_full.pdf".\
                            format(out_base,safe_name)
        Plotting.rupture_plot(true,pred,fig=fig)
        PlotUtilities.savefig(fig,final_out_rupture)
        fig = PlotUtilities.figure((7,4))
        # plot the metric plot
        Plotting.rupture_plot(true,pred,use_legend=True,
                              distance_histogram=distance_histogram,
                              fig=fig,color_pred=color_pred)
        final_out_path = "{:s}{:s}_full.pdf".format(out_base,safe_name)
        PlotUtilities.label_tom(fig,loc=locs)
        plt.suptitle(name,y=0.98,color=colors[i],alpha=0.7)
        PlotUtilities.save_twice(fig,final_out_path +".svg",
                                 final_out_path + ".png",
                                 subplots_adjust=dict(wspace=0.3,hspace=0.1,
                                                      left=0.1,bottom=0.1,
                                                      top=0.87))




if __name__ == "__main__":
    run()
