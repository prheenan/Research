# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from Research.Personal.EventDetection.OtherMethods import method_helper
from Research.Personal.EventDetection.Util import Analysis,Plotting,Scoring

from GeneralUtil.python import PlotUtilities
import fovea

def plot_classification(out_base,identifier,ex,scorer):
    """
    makes a plot of a classifier at (<out_base>+<identifier>)

    Args:
        out_base: directory to save
        identifier: name to save out
        ex: the split_fec we feed as the first arugment to funciton
        scorer: see single_classification

    Returns:
        Nothing
    """
    fig = PlotUtilities.figure()
    Plotting.plot_classification(ex,scorer)
    PlotUtilities.savefig(fig,out_base + identifier + ".png")

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    fec = method_helper.get_example()
    split_fec = Analysis.zero_and_split_force_extension_curve(fec)
    out_base = "./out/"
    scorers = []
    weights = np.linspace(start=0.1,stop=0.15,num=15)
    for w in weights:
        peaks_predicted = fovea.predict(fec,weight=w)
        scorer = Scoring.get_scoring_info(split_fec,peaks_predicted)
        plot_classification(out_base,"weight={:.3f}".format(w),split_fec,scorer)
        scorers.append(scorer)
    distances = [s.minimum_distance_median() for s in scorers]
    true = scorers[0].n_events
    predicted = np.array([s.n_events_predicted for s in scorers])
    relative_error = np.abs(predicted-true)/true
    xlim = [min(weights)*0.95,max(weights)*1.105]
    fig = PlotUtilities.figure(figsize=(8,12))
    n_plots = 2
    plt.subplot(n_plots,1,1)
    plt.plot(weights,distances,'bo-')
    plt.xlim(xlim)
    PlotUtilities.lazyLabel("","Median distance-to-events","")
    plt.subplot(n_plots,1,2)
    plt.plot(weights,relative_error,'bo-')
    plt.xlim(xlim)
    PlotUtilities.lazyLabel("Weights","Relative error in number of events","")
    PlotUtilities.savefig(fig,out_base+"weights.png")
        
if __name__ == "__main__":
    run()
