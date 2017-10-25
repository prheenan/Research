# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,cProfile,os,copy

sys.path.append("../../../../../../../../")
sys.path.append("../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from GeneralUtil.python import GenUtilities,CheckpointUtilities
import Util

def split(input_dir):
    """
    reads in the data, filters it, and decimates it 
    
    Args:
        input_dir: where the data to be read are
        n_filter_points: how many savitsky-golay points to use
        n_decimate: step size to use in the filtered
    Returns:
        yields the examples as it gets them 
    """
    # read in the data we want 
    examples = CheckpointUtilities.lazy_multi_load(input_dir)
    # filter (and decimate) the data as desired 
    for e in examples:
        approach,retract = FEC_Util.GetApproachRetract(e)
        yield retract 

def run():
    """
    Filters the input data to something manageable. 
    """
    input_dir = Util.cache_filtered()
    cache_dir = Util.cache_retract()
    GenUtilities.ensureDirExists(cache_dir)
    load_f = lambda: split(input_dir)
    CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_f,
                                   force=True,
                                   name_func=FEC_Util.fec_name_func)
    
 
        
if __name__ == "__main__":
    run()
