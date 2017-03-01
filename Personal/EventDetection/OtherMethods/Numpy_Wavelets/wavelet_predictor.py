# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from Research.Personal.EventDetection.Util import Analysis,Plotting,Scoring
from scipy.signal import find_peaks_cwt

def predict(fec,widths=None,min_snr=10,**kwargs):
    split_fec = Analysis.zero_and_split_force_extension_curve(fec)    
    force = split_fec.retract.Force
    if (widths is None):
        expected_max_log_width = int(np.ceil(np.log2(0.05 * force.size)))
        widths = np.logspace(0,expected_max_log_width,base=2,num=3)
    peak_indices = find_peaks_cwt(force,widths=widths,min_snr=min_snr,**kwargs)
    return peak_indices
