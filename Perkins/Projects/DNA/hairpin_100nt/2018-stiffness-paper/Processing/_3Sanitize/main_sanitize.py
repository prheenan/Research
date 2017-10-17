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

def sanitize(input_dir,seconds_for_zero):
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
        # get a copy of the data
        tmp = e._slice(slice(0,None,1))
        # determine where the zero force is 
        max_t = max(e.Time)
        idx_zero = np.where(e.Time < max_t - seconds_for_zero)[0]
        assert idx_zero.size > 0 
        # POST: something to look for 
        last_zero_idx = idx_zero[-1]
        e.Force -= np.median(e.Force[last_zero_idx])
        e.Force *= -1 
        # determine the zero crossing 
        zero_idx_sep = np.where(e.Force > 0)[0]
        assert zero_idx_sep.size > 0
        e.Separation -= e.Separation[zero_idx_sep[0]]
        yield e
        
def run():
    """
    Filters the input data to something manageable. 
    """
    seconds_for_zero = 0.25
    input_dir = Util.cache_retract()
    cache_dir = Util.cache_sanitized()
    GenUtilities.ensureDirExists(cache_dir)
    load_f = lambda: sanitize(input_dir,seconds_for_zero)
    CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_f,
                                   force=True,
                                   name_func=FEC_Util.fec_name_func)
    
 
        
if __name__ == "__main__":
    run()
