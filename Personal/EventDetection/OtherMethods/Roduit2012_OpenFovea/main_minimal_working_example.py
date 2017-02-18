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

def call_fovea(split_fec,weight=10):
    retract = split_fec.retract
    kwargs = dict(curve_x = retract.Separation*1e9,
                  curve_y = retract.Force*1e12,
                  # since we are already giving separation and force, we dont
                  # want to do any transformation
                  deflection_sensitivity = 1,
                  spring_constant = 1,
                  # assume that the point of contact is already figured out
                  poc=0,
                  weight=0.15,
                  # assume no drift, so baseline (coefficients for
                  # linear, drift-correcting fit)  is zero
                  baseline=[0,0],
                  fit_model=None)
    find = curve.event_find(**kwargs)
    mean_slice = lambda ev: int(np.round(np.mean([ev.start,ev.stop])))
    event_center = [mean_slice(f['Slice']) for f in find]
    return event_center

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
    time,separation,force = retract.Time,retract.Separation,retract.Force
    peaks_predicted = call_fovea(ex)
    # convert force to pN for this example
    scorer = Scoring.get_scoring_info(ex,peaks_predicted)
    fig = PlotUtilities.figure()
    Plotting.plot_classification(ex,scorer)
    PlotUtilities.savefig(fig,"./out_hooke.png")

if __name__ == "__main__":
    run()
