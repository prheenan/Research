# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle,os,string

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting,Analysis,\
    Learning
from Research.Personal.EventDetection.Util.Offline import plotting_metrics
import matplotlib.gridspec as gridspec

from Research.Personal.EventDetection.Util.Plotting \
    import algorithm_colors,algorithm_markers,algorithm_linestyles

def dist_values(data_file):
    trials = CheckpointUtilities.lazy_load(data_file)
    n_cols = 3
    to_ret = []
    for t in trials:
        ruptures_valid_true,ruptures_valid_pred = \
                Learning.get_true_and_predicted_ruptures_per_param(t)
        ret  = [Offline.metrics(true,pred) 
                for true,pred in zip(ruptures_valid_true,ruptures_valid_pred)]
        coeffs = [r[0] for r in ret]
        bcc = [1-c[-1] for c in coeffs]
        to_ret.append([bcc,t.param_values(),t.description.lower()])
    return to_ret

def run(base="./"):
    """
    
    """
    out_base = base
    data_file = "../FigurePerformance_CS/data/Scores.pkl"
    force=False
    events_off = CheckpointUtilities.getCheckpoint("cache_dist.pkl",
                                                   dist_values,force,
                                                   data_file)
    max_with_error = 0
    min_with_error = np.inf
    # get the bounds
    for i,(args,x,name) in enumerate(events_off):
        max_with_error = max(max_with_error,max(args))
        min_with_error = min(min_with_error,min(args))
    ylim = [min_with_error/2,max_with_error*2]
    # plot everything
    fig = PlotUtilities.figure((16,8))
    n_cols = 3
    parameters = ["FEATHER probability","OpenFovea sensitivity",
                  "Scientific Python minimum SNR"]
    markers = Plotting.algorithm_markers()
    colors = Plotting.algorithm_colors()
    for i,(args,x,name) in enumerate(events_off):
        first = (i == 0)
        param = parameters[i]
        xlabel = "Tuning Parameter\n({:s})".format(param)
        ylabel = "" if not first else "BCC"
        plt.subplot(1,n_cols,(i+1))
        lazy_kwargs = dict(useLegend=first,
                           frameon=True)
        plt.loglog(x,args,marker=markers[i],linestyle='-',color=colors[i],
                   markersize=7)
        plot_name = Plotting.algorithm_title_dict()[name]
        title = "Tuning curve for {:s}".format(plot_name)
        PlotUtilities.lazyLabel(xlabel,ylabel,title,**lazy_kwargs)
        PlotUtilities.tickAxisFont()
        plt.ylim([min_with_error/2,max_with_error*2])
    loc = (-0.10,1.025)
    PlotUtilities.label_tom(fig,loc=loc)
    PlotUtilities.savefig(fig,out_base + "tuning.pdf")


if __name__ == "__main__":
    run()
