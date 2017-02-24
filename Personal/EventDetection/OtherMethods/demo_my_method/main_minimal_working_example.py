# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from Research.Personal.EventDetection.OtherMethods import method_helper
from Research.Personal.EventDetection.Util import Analysis,Plotting,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector

from GeneralUtil.python import PlotUtilities

def adhesion_mask(surface_index,probability_distribution,threshold):
    """
    returns a boolean mask which is 0 where we can predict and adhesion 
    and zero elsewhere
    
    Args:
        surface_index: our best guess for where the surface is. 
        probability_distribution: see Detector._event_mask
        threshold: see Detector._event_mask
    Returns:
        list of event slices
    """
    to_ret = np.ones(probability_distribution.size,dtype=np.bool_)
    non_events = probability_distribution > threshold
    # remove all things before the predicted surface
    to_ret[:surface_index] = 0
    # remove everything until we arent at an event anymore
    non_events_after_predicted_surface = np.where(non_events[surface_index:])[0]
    if (non_events_after_predicted_surface.size == 0):
        # everything is always an event after the surface. Not much we can do,
        # so dont mess with to_ret
        pass
    else:
        # zero out everything until the first non-event
        first_non_event = non_events_after_predicted_surface[0]
        absolute_index_of_non_event = surface_index + first_non_event
        # XXX add in minimum distance between?
        # XXX add in 'need to not have an event for a certain amount of time'?
        to_ret[:absolute_index_of_non_event] = 0
    return to_ret
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    ex = method_helper.get_example()
    thresh = 0.02
    # XXX : get first index where less than the threshold after surface index?
    surface_index = ex.get_predicted_retract_surface_index()
    adhesion_func = lambda *args: adhesion_mask(surface_index,*args)
    info = Detector._predict_helper(ex,threshold=thresh,
                                    condition_function=adhesion_func)
    # XXX fix threshhold
    fig = PlotUtilities.figure(figsize=(8,12))    
    Plotting.plot_prediction_info(ex,info)
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
