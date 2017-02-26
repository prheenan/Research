# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from Research.Personal.EventDetection.Util import InputOutput,Analysis,Plotting
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis \
    import FEC_Util

def get_example(base="./"):
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
    file_name = base +\
                "2017-02-04-masters-data-650nm-dna-2.5ng_ul_1-15_dilution-25"+\
                "-hours-in-pbs-multiple-loading-rates-170203-500nm-s.pxp_" + \
                "Image1849Concat.csv"
    file_path = file_folder_path + file_name
    example = InputOutput.read_and_cache_file(file_path,cache_directory=base,
                                              force=False,has_events=True)
    example_split = Analysis.zero_and_split_force_extension_curve(example)
    return example_split
