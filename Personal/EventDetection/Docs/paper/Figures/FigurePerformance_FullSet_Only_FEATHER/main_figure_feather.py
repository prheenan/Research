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
    out_base = in_base
    data_file = in_base + "data/Scores.pkl"
    force=False
    cache_file = out_base + "cache.pkl"
    metrics = CheckpointUtilities.getCheckpoint(cache_file,get_best_metrics,
                                                False,data_file)
    loc_left = (-0.10,1.1)
    loc_top = (-0.12,1.05)
    loc_lower = (-0.12,0.95)
    locs = [loc_top,loc_left,loc_left,loc_lower,loc_lower]
    titles = Plotting.algorithm_title_dict()
    colors = Plotting.algorithm_colors()
    for i,m in enumerate(metrics):
        name = titles[m.name.lower()]
        fig = PlotUtilities.figure((16,8))
        distance_histogram= Offline.event_error_kwargs(m)
        true,pred = m.true,m.pred
        # plot the metric plot
        Plotting.rupture_plot(true,pred,use_legend=True,
                              distance_histogram=distance_histogram,
                              fig=fig)
        final_out_path = "{:s}{:s}.pdf".format(out_base,name)
        PlotUtilities.label_tom(fig,loc=locs,fontsize=18)
        plt.suptitle(name,fontsize=25,y=0.95,color=colors[i],alpha=0.7)
        PlotUtilities.savefig(fig,final_out_path,
                              subplots_adjust=dict(wspace=0.2,hspace=0.1,
                                                   left=0.05,top=0.85))




if __name__ == "__main__":
    run()
