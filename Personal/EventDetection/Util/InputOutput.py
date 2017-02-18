# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os
from scipy import interpolate
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from scipy.stats import norm
from GeneralUtil.python import CheckpointUtilities


def read_and_cache_file(file_path,cache_directory,has_events=False,force=True):
    """
    reads the (csv) file at file_path, cachine it to cache_directory,
    reading in events 
    
    Args;
        file_path: where the file lives
        cache_directory: where to cache the file
        has_events:  if this file has events
        force: if true, force a read
    Returns
        TimeSepForce Object
    """
    file_name = os.path.basename(file_path)
    cache_file = cache_directory + file_name+ ".pkl"
    func_to_call = FEC_Util.read_time_sep_force_from_csv
    return CheckpointUtilities.getCheckpoint(cache_file,func_to_call,force,
                                             file_path,has_events=has_events)
                   