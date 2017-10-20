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
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
from IgorUtil.PythonAdapter import TimeSepForceObj
import Util
from FitUtil.FreelyJointedChain.Python.Code import FJC

    
def fit_polymer_model(example):
    """
    Fits the polymer model de jour to example
    
    Args:
        example: the timesepforce to fit to
    Returns:
        x0, the parameters of the fit 
    """
    wlc_params = dict(K0=2000e-12,kbT=4.1e-21,Lp=0.3e-9)
    ranges = [(10e-9,100e-9)]
    fit_dict = dict(brute_dict=dict(Ns=50,ranges=ranges),
                    **wlc_params)
    x_raw,y_raw = example.Separation,example.Force
    x0,model_x,model_y = FJC.fit_fjc_contour(x_raw,y_raw,**fit_dict)        
    return x0,model_x,model_y


def align(input_dir,m_to_fit_before_max=20e-9,arbitrary_offset_m=-90e-9):
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
        # determine where the max force is
        max_idx = np.argmax(tmp.Force)
        sep_max = tmp.Separation[max_idx]
        idx_where_le = \
            np.where(tmp.Separation < sep_max - m_to_fit_before_max)[0]
        assert idx_where_le.size > 0 
        idx_fit = idx_where_le[-1]
        obj_to_fit = tmp._slice(slice(idx_fit,max_idx))
        x0,model_x,model_y = fit_polymer_model(obj_to_fit)
        L0 = x0
        tmp.Separation -= (L0+arbitrary_offset_m)
        yield tmp
    
        
def run():
    """
    Filters the input data to something manageable. 
    """
    m_to_fit_before_max = 20e-9
    input_dir = Util.cache_sanitized()
    cache_dir = Util.cache_aligned()
    GenUtilities.ensureDirExists(cache_dir)
    load_f = lambda: align(input_dir,m_to_fit_before_max)
    data = CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_f,
                                          force=False,
                                          name_func=FEC_Util.fec_name_func)
    fig = PlotUtilities.figure()
    plt.subplot(1,1,1)
    FEC_Plot.heat_map_fec(data,num_bins=(100,500),separation_max=100)
    plt.ylim([-10,None])
    PlotUtilities.savefig(fig,"./out.png")
 
        
if __name__ == "__main__":
    run()
