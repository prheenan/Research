# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../")
sys.path.append("../../../../../../../../../")
from LandscapeGeneration import LandscapeUtil
from Processing import Util
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot

from scipy.interpolate import LSQUnivariateSpline


def get_aligned_regions(input_dir):
    """
    Returns the 'aligned' iwt regions 
    """
    examples = CheckpointUtilities.lazy_multi_load(input_dir)
    spline_t = [e.spline(e.data.Time) for e in examples]
    # get the bounds for the region
    max_of_min = max([min(s) for s in spline_t])
    # XXX why is this not min(max( ? 
    min_of_max = min([max(s) for s in spline_t])    
    for i,example_info in enumerate(examples):
        e = example_info.data
        z = spline_t[i]
        where_between_idx = np.where( (z <= max_of_min) &
                                      (z >= min_of_max))[0]
        slice_idx = example_info._idx_pairs                                      
        # determine the regions we care about 
        t = e.Time
        plt.plot(t,z)
        plt.axhline(max_of_min)
        plt.axhline(min_of_max)
        for s in slice_idx:
            plt.axvline(t[s.start])
        plt.show()

def run():
    input_dir =  LandscapeUtil.cache_landscape_regions("../../../")
    cache_dir = LandscapeUtil.cache_landscape_aligned_regions("../../../")
    load_f = lambda: get_aligned_regions(input_dir)
    e = CheckpointUtilities.multi_load(cache_dir=cache_dir,
                                       load_func=load_f,
                                       force=True,limit=None,
                                       name_func=FEC_Util.fec_name_func)

if __name__ == "__main__":
    run()
