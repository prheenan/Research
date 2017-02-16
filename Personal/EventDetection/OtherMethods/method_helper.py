# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from Research.Personal.EventDetection.Util import InputOutput,Analysis,Plotting
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis \
    import FEC_Util

def get_example():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    file_folder_path = FEC_Util.default_data_root() + \
        "4Patrick/CuratedData/Masters_CSCI/Positive/650nm-4x-bio/csv/" + \
        "500-nanometers-per-second/"
    file_name = "2017-02-04-masters-data-650nm-dna-2.5ng_ul_1-15_dilution-25"+\
                "-hours-in-pbs-multiple-loading-rates-170203-500nm-s.pxp_" + \
                "Image1849Concat.csv"
    file_path = file_folder_path + file_name
    example = InputOutput.read_and_cache_file(file_path,cache_directory="./",
                                              force=False,has_events=True)
    example_split = Analysis.split_FEC_by_meta(example)
    approach = example_split.approach
    retract = example_split.retract 
    # get the autocorrelation time of the retract force (what we care about)
    x,f = retract.Time,retract.Force
    dx = np.median(np.diff(x))
    frequency = retract.ThermalFrequency
    n_points_fit = int(np.ceil((1/frequency) * (1/dx)))
    tau,auto_coeffs,auto_correlation = \
        Analysis.auto_correlation_tau(x,f,deg_autocorrelation=1,
                                      fit_idx_max=n_points_fit)
    num_points = int(np.ceil(tau/dx))
    # zero out everything to the approach using the autocorrelation time 
    Analysis.zero_by_approach(example_split,num_points)                                              
    return example_split