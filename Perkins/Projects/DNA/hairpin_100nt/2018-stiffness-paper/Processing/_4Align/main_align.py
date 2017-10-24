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
from scipy.interpolate import LSQUnivariateSpline

    
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


def align(input_dir,m_to_fit=20e-9):
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
        t = tmp.Time
        force = tmp.Force
        n_bins = int(np.ceil(force.size/100))
        knots = np.linspace(min(t),max(t),endpoint=True,num=n_bins)
        spline_force = LSQUnivariateSpline(x=t,y=tmp.Force,t=knots[1:-1],
                                           k=3)
        spline_sep = LSQUnivariateSpline(x=t,y=tmp.Separation,t=knots[1:-1],
                                         k=3)
        spline_force_t = spline_force(t)
        spline_sep_t = spline_sep(t)
        where_ge_0 = np.where(spline_force_t > 0)[0]
        assert where_ge_0.size > 0 , "Data not zeroed."
        zero_idx = where_ge_0[0]
        max_force_idx = np.argmax(spline_force_t)
        sep_until_max_force = spline_sep_t[:max_force_idx]
        idx_where_le = np.where(sep_until_max_force <= \
                                max(sep_until_max_force) - m_to_fit)[0]
        assert idx_where_le.size > 0 
        idx_fit_i = idx_where_le[-1]
        idx_fit_f = max_force_idx
        slice_fit = slice(idx_fit_i,idx_fit_f)
        obj_to_fit = tmp._slice(slice_fit)
        zero_sep = spline_sep_t[zero_idx]
        obj_to_fit.Separation -= zero_sep
        x0,model_x,model_y = fit_polymer_model(obj_to_fit)
        L0 = x0
        min_sep,min_z = min(tmp.Separation),min(tmp.ZSnsr)
        offset_sep = 0
        offset_z = 0
        to_ret = e._slice(slice(zero_idx,None,1))
        to_ret.Separation -= offset_sep
        to_ret.ZSnsr -= offset_z
        yield to_ret
    
        
def run():
    """
    Filters the input data to something manageable. 
    """
    m_to_fit = 30e-9
    input_dir = Util.cache_sanitized()
    cache_dir = Util.cache_aligned()
    GenUtilities.ensureDirExists(cache_dir)
    load_f = lambda: align(input_dir,m_to_fit)
    data = CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_f,
                                          force=True,
                                          name_func=FEC_Util.fec_name_func)
    for d in data:
        plt.subplot(2,1,1)
        plt.plot(d.Separation,d.Force)
        plt.subplot(2,1,2)
        plt.plot(d.Time,d.Force)
    plt.show()
    fig = PlotUtilities.figure()
    plt.subplot(1,1,1)
    FEC_Plot.heat_map_fec(data,num_bins=(100,500),separation_max=100)
    plt.ylim([-10,None])
    PlotUtilities.savefig(fig,"./out.png")
 
        
if __name__ == "__main__":
    run()
