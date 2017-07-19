# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from Research.Personal.EventDetection.OtherMethods.Roduit2012_OpenFovea \
    import curve
    
    
from Research.Personal.EventDetection.Util import Analysis

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
         location of event centers as a list.
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

def predict(fec,weight=0.1):
    """
    Calls the open fovea method of event detection 

    Args:
         fec: the TimeSepForce object we want to use

         weight: see call_fovea
    Returns:
         location of event centers as a list.
    """
    split_fec = Analysis.zero_and_split_force_extension_curve(fec)
    surface_index_retract = split_fec.get_predicted_retract_surface_index()
    return call_fovea(split_fec,weight=weight,poc=surface_index_retract)
