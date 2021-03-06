# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from Research.Personal.EventDetection.Util.Learning import _get_single_curve
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from Research.Personal.EventDetection.OtherMethods.Roduit2012_OpenFovea import \
    fovea
from Research.Personal.EventDetection.OtherMethods.Numpy_Wavelets import \
    wavelet_predictor
    

def get_learners(n_points_no_event=5,n_points_fovea=5,n_points_wavelet=5,
                 no_event_log10_start=-3,no_event_log10_end=-2):
    """
    Returns a list of learning_curve objects

    Args:
        n_points_no_event: number of points for varying the no event portion of 
        things
            
        n_points_fovea: number of points to use on fovea
    Returns:
        list of learning curves
    """
    # make the no event example
    no_event_func = lambda arg_list: [dict(threshold=t) for t in arg_list]
    no_event_tuple = [Detector.predict,np.logspace(no_event_log10_start,
                                                   no_event_log10_end,
                                                   endpoint=True,
                                                   base=10,
                                                   num=n_points_no_event)]
    no_event_curve = _get_single_curve("No Event",no_event_tuple,no_event_func) 
    # make the fovea example
    fovea_func = lambda arg_list: [dict(weight=w) for w in arg_list]
    fovea_tuple = [fovea.predict,np.logspace(-2,np.log10(0.5),
                                             endpoint=True,
                                             num=n_points_fovea)]
    fovea_curve = _get_single_curve("Open Fovea",fovea_tuple,fovea_func)
    # make the CWT example
    cwt_func = lambda arg_list: [dict(min_snr=w) for w in arg_list]
    cwt_tuple = [wavelet_predictor.predict,
                 np.logspace(start=np.log10(10),stop=np.log10(200),
                             num=n_points_wavelet)]
    wavelet_curve = _get_single_curve("Wavelet transform",cwt_tuple,cwt_func)   
    return [no_event_curve,fovea_curve,wavelet_curve]
