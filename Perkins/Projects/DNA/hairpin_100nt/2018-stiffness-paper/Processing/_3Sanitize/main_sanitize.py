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
from scipy.interpolate import LSQUnivariateSpline

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
        max_t = max(tmp.Time)
        idx_zero = np.where(tmp.Time < max_t - seconds_for_zero)[0]
        assert idx_zero.size > 0 
        # POST: something to look for 
        last_zero_idx = idx_zero[-1]
        tmp.Force -= np.median(tmp.Force[last_zero_idx])
        tmp.Force *= -1 
        t = e.Time
        n_bins = 200
        time_interp = np.linspace(min(t),max(t),endpoint=True,num=n_bins)
        spline_kw = dict(x=t,t=time_interp[1:-1],k=3)
        force_spline_fit = LSQUnivariateSpline(y=tmp.Force,**spline_kw)
        force_spline_t = force_spline_fit(t)
        where_ge_0 = np.where(force_spline_t > 0)[0]
        assert where_ge_0.size > 0 
        first_ge_0 = where_ge_0[0]
        sep_spline = LSQUnivariateSpline(y=e.Separation,**spline_kw)(t)
        z_spline = LSQUnivariateSpline(y=e.ZSnsr,**spline_kw)(t)
        tmp.Separation -= sep_spline[first_ge_0]
        tmp.ZSnsr -= z_spline[first_ge_0]
        yield tmp
        
def run():
    """
    Filters the input data to something manageable. 
    """
    seconds_for_zero = 0.25
    input_dir = Util.cache_retract()
    cache_dir = Util.cache_sanitized()
    GenUtilities.ensureDirExists(cache_dir)
    load_f = lambda: sanitize(input_dir,seconds_for_zero)
    e = CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_f,
                                       force=True,
                                       name_func=FEC_Util.fec_name_func)
    for tmp in e:
        plt.plot(tmp.Separation,tmp.Force)
    plt.show()
    
 
        
if __name__ == "__main__":
    run()
