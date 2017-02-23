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
from Research.Personal.EventDetection.OtherMethods.Roduit2012_OpenFovea.\
    openfovea_src.openfovea.fovea_toolbox import curve

def call_fovea(split_fec,weight=10,poc=0):
    """
    Calls the open fovea method of event detection 

    Args:
         split_fec: the Analysis.split_force_extension object to event find
         weight: see openfovea.fovea_toolbox.curve.event_find. Essentially, 
         a sensitivity parameter

         poc: see openfovea.fovea_toolbox.curve.event_find. The index in the 
         retract where the surface is.
    Returns:
         locaiton of event centers as a list.
    """
    retract = split_fec.retract
    kwargs = dict(curve_x = retract.Separation*1e9,
                  curve_y = retract.Force*1e12,
                  # since we are already giving separation and force, we dont
                  # want to do any transformation
                  deflection_sensitivity = 1,
                  spring_constant = 1,
                  # poc: point of contact. 
                  poc=poc,
                  weight=weight,
                  # assume no drift, so baseline (coefficients for
                  # linear, drift-correcting fit)  is zero
                  baseline=[0,0],
                  fit_model=None)
    find = curve.event_find(**kwargs)
    mean_slice = lambda ev: int(np.round(np.mean([ev.start,ev.stop])))
    if (find is not None):
        event_center = [mean_slice(f['Slice']) for f in find]
    else:
        event_center = []
    return event_center

def single_classification(ex,function,**kwargs):
    """
    calls a single classification function, getting the scorer

    Args:
        ex: the split_fec we feed as the first arugment to funciton
        function: function taking in ex as the first argument, then other args,
        returning a list of event centers.
        
        **kwargs: passed to function after ex
    Returns:
        Scoring.score object for this classifier on this data
    """
    peaks_predicted = function(ex,**kwargs)
    # convert force to pN for this example
    scorer = Scoring.get_scoring_info(ex,peaks_predicted)
    return scorer

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
    ex = method_helper.get_example()
    retract = ex.retract
    _,surface_index,_ = Analysis.get_surface_index(retract,
                                                   n_smooth=ex.tau_num_points,
                                                   last_less_than=False)
    time,separation,force = retract.Time,retract.Separation,retract.Force
    out_base = "./out/"
    scorers = []
    weights = np.linspace(start=0.1,stop=0.15,num=15)
    for w in weights:
        kwargs_fovea = dict(weight=w,
                            poc=surface_index)
        scorer = single_classification(ex,call_fovea,**kwargs_fovea)
        plot_classification(out_base,"weight={:.3f}".format(w),ex,scorer)
        scorers.append(scorer)
    distances = [s.minimum_distance_median for s in scorers]
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
