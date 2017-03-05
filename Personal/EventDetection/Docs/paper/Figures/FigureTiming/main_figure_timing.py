# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import Learning,InputOutput
from Research.Personal.EventDetection.Timing import TimePlot
from Research.Personal.EventDetection.Timing.Timer import \
    time_trials_by_loading_rate,time_trials


def run(base="./"):
    """
    
    """
    # XXX should be able to get curve numbers of out pickle
    curve_numbers = [1,2,5,10,30,50,100,150,200]
    data_file = base + "all.pkl"
    with open(data_file,"r") as f:
        times = pickle.load(f)
    out_base = base + "out/"
    GenUtilities.ensureDirExists(out_base)
    # sort the times by their loading rates
    max_time = max([l.max_time_trial() for l in times])
    min_time = min([l.min_time_trial() for l in times])
    # plot the Theta(n) coefficient for each
    fig = PlotUtilities.figure()
    TimePlot.plot_learner_prediction_time_comparison(times)
    PlotUtilities.legend(loc="lower right",frameon=True)
    PlotUtilities.savefig(fig,out_base + "compare.png")
    for learner_trials in times:
        base_name = out_base + learner_trials.learner.description
        # plot the timing veruses loading rate and number of points 
        fig = PlotUtilities.figure()
        TimePlot.plot_learner_versus_loading_rate_and_number(learner_trials)
        fudge_x_low = 10
        fudge_x_high = 2
        fudge_y = 1.5
        plt.ylim([min_time/fudge_y,max_time*fudge_y])
        plt.xlim([1/fudge_x_low,max(curve_numbers)*fudge_x_high])
        plt.yscale('log')
        plt.xscale('log')        
        PlotUtilities.legend(loc="upper left",frameon=True)
        PlotUtilities.savefig(fig,  base_name + "_all_trials.png")
        # plot the slopes
        fig = PlotUtilities.figure()
        TimePlot.plot_learner_slope_versus_loading_rate(learner_trials)
        PlotUtilities.legend(loc="lower right",frameon=True)
        PlotUtilities.savefig(fig, base_name + "_slopes.png")


if __name__ == "__main__":
    run()
