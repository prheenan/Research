# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,cProfile,os

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from FitUtil.FreelyJointedChain.Python.Code import FJC
from GeneralUtil.python import PlotUtilities


def hairpin_plots(example,filter_fraction,out_path):
    n_filter = int(np.ceil(example.Force.size * filter_fraction))
    # make a plot vs time 
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    FEC_Plot.force_versus_time(example,NFilterPoints=n_filter)
    plt.show()
    plt.subplot(2,1,2)
    FEC_Plot.z_sensor_vs_time(example,NFilterPoints=n_filter)
    PlotUtilities.savefig(fig,out_path + "vs_time.png")
    # make a force-extension plot
    fig = PlotUtilities.figure()
    FEC_Plot.FEC(example,NFilterPoints=n_filter)
    PlotUtilities.savefig(fig,out_path + "vs_sep.png")
        
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    abs_dir = "./"
    examples = FEC_Util.\
        cache_individual_waves_in_directory(pxp_dir=abs_dir,force=False,
                                            cache_dir="./cache/",limit=1)
    ### XXX TODO
    # (1) correct for interference artifact
    # (2) get regions for WLC fit
    # (3) fit WLC to regions
    # (4) Invert WLC, determine dsDNA and ssDNA contour lengths at each force 
    region_fit = [8.1,8.9]
    
    for i,ex in enumerate(examples):
        _,processed = FEC_Util.PreProcessFEC(ex)
        # XXX be more clever.
        sep_offset = 13e-9
        force_offset = 8e-12
        processed.Separation -= sep_offset
        processed.Force -= force_offset
        for_wlc_fit = FEC_Util.slice_by_time(processed,*region_fit) 
        wlc_params = dict(K0=2000e-12,Lp=0.2e-9,kbT=4.1e-21)
        fit_dict = dict(brute_dict=dict(Ns=20,ranges=[(10e-9,90e-9)]),
                        **wlc_params)
        x_raw,y_raw = for_wlc_fit.Separation,for_wlc_fit.Force
        x0,model_x,model_y = FJC.fit_fjc_contour(x_raw,y_raw,**fit_dict)
        print(x0)
        plt.plot(x_raw,y_raw)
        plt.plot(model_x,model_y)
        plt.show()
    return 
    for i,ex in enumerate(examples):
        hairpin_plots(ex,filter_fraction=1e-3,out_path="./out/{:d}".format(i))
        
if __name__ == "__main__":
    run()
