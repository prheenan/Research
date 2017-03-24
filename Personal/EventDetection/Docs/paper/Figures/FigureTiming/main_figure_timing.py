# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import Learning,InputOutput,Plotting
from Research.Personal.EventDetection.Timing import TimePlot
from Research.Personal.EventDetection.Timing.Timer import \
    time_trials_by_loading_rate,time_trials

from Research.Personal.EventDetection.Util.Plotting \
    import algorithm_colors,algorithm_markers,algorithm_linestyles


def get_supplemental_figure(output_path,trials):
    """
    creates the 'supplemental' timing figure (for use in the appendix) 

    Args:
        see get_main_figure
    """
    # XXX should be able to get curve numbers of out pickle
    curve_numbers = [1,2,5,10,30,50,100,150,200]
    # make a plot comparing *all* of the Big-O plots of the data
    plt.subplot(3,2,1)
    # sort the times by their loading rates
    max_time = max([l.max_time_trial() for l in trials])
    min_time = min([l.min_time_trial() for l in trials])
    fig = PlotUtilities.figure((16,16))
    # plot the Theta(n) coefficient for each
    n_rows = 3
    n_cols = 2
    style_data,style_pred = style_data_and_pred()
    x_label = "C (number of curves)"
    y_label = "Runtime (s)"
    x_label_big_o = "N (points per curve)"
    y_label_big_o = "Runtime per curve (s) "
    ylim_big_o = [1e-3,1e3]
    for i,learner_trials in enumerate(trials):
        description = TimePlot.learner_name(learner_trials)
        plot_idx = i*n_cols+1
        plt.subplot(n_rows,n_cols,plot_idx)
        # plot the timing veruses loading rate and number of points 
        TimePlot.plot_learner_versus_loading_rate_and_number(learner_trials)
        fudge_x_low = 20
        fudge_x_high = 2
        fudge_y = 4
        plt.ylim(ylim_big_o)
        plt.xlim([1/fudge_x_low,max(curve_numbers)*fudge_x_high])
        plt.yscale('log')
        plt.xscale('log')        
        useLegend= (i == 0)
        last = (i == (len(trials) - 1))
        PlotUtilities.lazyLabel("","","",useLegend=useLegend,frameon=True,
                                legend_kwargs=dict(fontsize=15))
        if (not useLegend):
            plt.gca().legend().remove()
        PlotUtilities.ylabel(y_label)
        if (last):
            PlotUtilities.xlabel(x_label)            
        PlotUtilities.title("Total runtime ({:s})".\
                            format(description))
        plt.subplot(n_rows,n_cols,plot_idx+1)
        style_dict = dict(style_data=style_data[i],style_pred=style_pred[i])
        TimePlot.plot_learner_slope_versus_loading_rate(learner_trials,
                                                        **style_dict)
        PlotUtilities.title("Runtime/curve of length N ({:s})".\
                            format(description))
        if (last):
            PlotUtilities.xlabel(x_label_big_o)
        else:
            PlotUtilities.xlabel("")
        PlotUtilities.ylabel(y_label_big_o)
        plt.ylim(ylim_big_o)
    PlotUtilities.label_tom(fig,loc=(-0.15,1.05))
    PlotUtilities.savefig(fig, output_path)


    # make a plot comparing the constants
    pass

def style_data_and_pred():
    colors,markers,linestyles = algorithm_colors(),algorithm_markers(),\
                                algorithm_linestyles()
    style_data=[dict(marker=m,color=c,linestyle='None')
                for m,c in zip(markers,colors)]
    style_pred = []
    for i in range(len(colors)):
        style_pred.append(dict(color=colors[i],linestyle=linestyles[i]))
    return style_data,style_pred
    

def make_main_figure(output_path,trials):
    """
    creates the 'main' timing figure (for use in the paper) 

    Args:
        output_path: where to save the file
        trials: the pickled timing trials information
    """
    xlim = [500,1e6]
    ylim=[1e-2,5e2]
    style_data,style_pred = style_data_and_pred()
    colors = algorithm_colors()
    # picturing 3x2, where we show the 'in the weeds' plots...
    fig = PlotUtilities.figure(figsize=(16,6))
    plt.subplot(1,2,1)
    for i,learner_trials in enumerate(trials):
        TimePlot.\
            plot_learner_slope_versus_loading_rate(learner_trials,
                                                   style_data=style_data[i],
                                                   style_pred=style_pred[i])
    PlotUtilities.legend(loc="upper left",frameon=True)
    PlotUtilities.title("")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.subplot(1,2,2)
    TimePlot.plot_learner_prediction_time_comparison(trials,color=colors)
    plt.ylim([1e3,3e6])
    PlotUtilities.legend(loc="lower right",frameon=True)
    PlotUtilities.title("")
    PlotUtilities.label_tom(fig,loc=(-0.05,1))
    PlotUtilities.savefig(fig,output_path)

def run(base="./"):
    """
    
    """
    data_file = base + "all.pkl"
    with open(data_file,"r") as f:
        trials = pickle.load(f)
    out_base = base 
    GenUtilities.ensureDirExists(out_base)
    make_main_figure(output_path=out_base+"timing.pdf",trials=trials)
    get_supplemental_figure(output_path=out_base+"supplemental.pdf",
                            trials=trials)

if __name__ == "__main__":
    run()
