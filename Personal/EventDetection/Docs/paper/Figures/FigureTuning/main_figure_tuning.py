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

def dist_values(data_file):
    trials = CheckpointUtilities.lazy_load(data_file)
    n_cols = 3
    to_ret = []
    for i,m in enumerate(trials):
        x_values = m.param_values()
        train_scores = m._scores_by_params(train=True)
        valid_scores = m._scores_by_params(train=False)
        tmp = Plotting.get_train_test_n_off_and_error(x_values,
                                                      train_scores,
                                                      valid_scores)
        to_ret.append([tmp,m.description.lower()])
    return to_ret
        

def run(base="./"):
    """
    
    """
    out_base = base
    data_file = base + "data/Scores.pkl"
    force=False
    events_off = CheckpointUtilities.getCheckpoint("cache_dist.pkl",
                                                   dist_values,force,
                                                   data_file)
    max_with_error = 0
    min_with_error = np.inf
    # get the bounds
    for i,(args,name) in enumerate(events_off):
        x_train,train_dist,train_dist_std,x_valid,valid_dist,valid_dist_std = \
            args
        local_max = max(max(train_dist+train_dist_std),
                        max(valid_dist+valid_dist_std))
        local_min = min(min(train_dist-train_dist_std),
                        min(valid_dist-valid_dist_std))
        max_with_error = max(max_with_error,local_max)
        min_with_error = min(min_with_error,local_min)
    ylim = [min_with_error/2,max_with_error*2]
    # plot everything
    fig = PlotUtilities.figure((16,8))
    n_cols = 3
    parameters = ["FEATHER probability","OpenFovea sensitivity",
                  "Scientific Python minimum SNR"]
    xlim = [None,[0.5,2],[10,200]]
    for i,(args,name) in enumerate(events_off):
        first = (i == 0)
        param = parameters[i]
        xlabel = "Tuning Parameter\n({:s})".format(param)
        ylabel = "" if not first else \
                 "Relative number of extra or missing events"
        plt.subplot(1,n_cols,(i+1))
        lazy_kwargs = dict(useLegend=first,
                           frameon=True)
        Plotting._plot_num_events_off(*args,xlabel=xlabel,ylabel=ylabel,
                                      lazy_kwargs=lazy_kwargs)
        plt.ylim(ylim)
        if (xlim[i] is not None):
            plt.xlim(xlim[i])
        PlotUtilities.tickAxisFont()
        plot_name = Plotting.algorithm_title_dict()[name]
        PlotUtilities.title("Tuning curve for {:s}".format(plot_name))
    loc = (-0.20,1.025)
    PlotUtilities.label_tom(fig,loc=loc)
    PlotUtilities.savefig(fig,out_base + "tuning.pdf")


if __name__ == "__main__":
    run()
