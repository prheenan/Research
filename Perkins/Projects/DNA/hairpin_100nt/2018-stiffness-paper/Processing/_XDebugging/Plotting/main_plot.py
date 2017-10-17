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

sys.path.append("../../../../../../../../../")
sys.path.append("../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
import Util

def _check_data(*data):
    key = data[0]
    for d in data:
        assert len(key) == len(d) , "Data not all the same length of date"
        # POST: same number
        for key_tmp,key_d in zip(key,d):
            assert FEC_Util.fec_name_func(0,key_tmp) == \
                   FEC_Util.fec_name_func(0,key_d) , "Name mismatch"
            # POST: names match
    # POST: all sets match 

def run():
    """
    Filters the input data to something manageable. 
    """
    dir_raw = Util.cache_raw("../../../")
    dir_filtered = Util.cache_filtered("../../../")
    out_dir = "./out/"
    NFilterPoints = 500
    GenUtilities.ensureDirExists(out_dir)
    # sort by the names...
    m_sort = lambda x : sorted(x,key=lambda e: FEC_Util.fec_name_func(0,e))
    examples_unfiltered = m_sort(CheckpointUtilities.lazy_multi_load(dir_raw))
    examples_filtered = \
        m_sort(CheckpointUtilities.lazy_multi_load(dir_filtered))
    _check_data(examples_unfiltered,examples_filtered)
    for ex_raw,ex_filt in zip(examples_unfiltered,examples_filtered):
        fig = PlotUtilities.figure()
        plt.subplot(3,1,1)
        plt.plot(ex_raw.Time,ex_raw.Force,color='k',alpha=0.3,label="raw")
        plt.plot(ex_filt.Time,ex_filt.Force,color='r',label="filtered")
        PlotUtilities.lazyLabel("Time (s)","Force (N)","")
        plt.subplot(3,1,2)
        FEC_Plot.FEC(ex_filt,PreProcess=True,NFilterPoints=1,
                     PlotLabelOpts=dict(linewidth=0.5))   
        PlotUtilities.lazyLabel("","Force (pN)","")
        plt.subplot(3,1,3)                        
        FEC_Plot.FEC(ex_raw,PreProcess=True,NFilterPoints=NFilterPoints,
                     PlotLabelOpts=dict(linewidth=0.5))  
        PlotUtilities.lazyLabel("Time (s)","Force (N)","")
        out_name = out_dir + FEC_Util.fec_name_func(0,ex_filt) + ".jpeg"
        PlotUtilities.savefig(fig,out_name)
    
 
        
if __name__ == "__main__":
    run()
