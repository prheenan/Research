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

def get_supplemental_figure(output_path,trials):
    """
    creates the 'supplemental' timing figure (for use in the appendix) 

    Args:
        see get_main_figure
    """

    # make a plot comparing *all* of the Big-O plots of the data

    # make a plot comparing the constants
    pass

def get_main_figure(output_path,trials):
    """
    creates the 'main' timing figure (for use in the paper) 

    Args:
        output_path: where to save the file
        trials: the pickled timing trials information
    """
    xlim = [500,1e6]
    ylim=[1e-2,5e2]
    # picturing 3x2, where we show the 'in the weeds' plots...
    colors = ['r','b','k']
    predict_styles = [dict(linestyle='--'),dict(linestyle='-.'),
                      dict(linestyle='-')]
    kw = dict(linestyle='None')
    style_data=[dict(marker='s',**kw),dict(marker='o',**kw),
                dict(marker='v',**kw)]
    markers = [s['marker'] for s in style_data]
    for i,c in enumerate(colors):
        style_data[i]['color'] = c
    fig = PlotUtilities.figure(figsize=(16,8))
    plt.subplot(1,2,1)
    for i,learner_trials in enumerate(trials):
        style_pred = dict(color=colors[i],**predict_styles[i])
        TimePlot.\
            plot_learner_slope_versus_loading_rate(learner_trials,
                                                   style_data=style_data[i],
                                                   style_pred=style_pred)
    PlotUtilities.legend(loc="upper left",frameon=True)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.subplot(1,2,2)
    TimePlot.plot_learner_prediction_time_comparison(trials,color=colors)
    PlotUtilities.legend(loc="lower right",frameon=True)
    PlotUtilities.title(r"No event asymptotic T(N) $\geq$10x faster")
    PlotUtilities.savefig(fig,output_path)


def run(base="./"):
    """
    
    """
    # XXX should be able to get curve numbers of out pickle
    curve_numbers = [1,2,5,10,30,50,100,150,200]
    data_file = base + "all.pkl"
    with open(data_file,"r") as f:
        trials = pickle.load(f)
    out_base = base + "out/"
    GenUtilities.ensureDirExists(out_base)
    get_main_figure(output_path=out_base+"timing.pdf",trials=trials)
    exit(1)
    # sort the times by their loading rates
    max_time = max([l.max_time_trial() for l in times])
    min_time = min([l.min_time_trial() for l in times])
    # plot the Theta(n) coefficient for each
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



if __name__ == "__main__":
    run()
