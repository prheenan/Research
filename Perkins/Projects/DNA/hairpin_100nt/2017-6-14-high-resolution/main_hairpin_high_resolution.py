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
from Research.Personal.EventDetection.Util import Analysis


def hairpin_plots(example,filter_fraction,out_path):
    n_filter = int(np.ceil(example.Force.size * filter_fraction))
    # make a plot vs time 
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    FEC_Plot.force_versus_time(example,NFilterPoints=n_filter)
    plt.subplot(2,1,2)
    FEC_Plot.z_sensor_vs_time(example,NFilterPoints=n_filter)
    PlotUtilities.savefig(fig,out_path + "vs_time.png")
    # make a force-extension plot
    fig = PlotUtilities.figure()
    FEC_Plot.FEC(example,NFilterPoints=n_filter)
    PlotUtilities.savefig(fig,out_path + "vs_sep.png")
        
def fit_polymer_model(example):
    wlc_params = dict(K0=2000e-12,kbT=4.1e-21)
    ranges = [(10e-9,90e-9),(0.1e-9,1e-9)]
    fit_dict = dict(brute_dict=dict(Ns=20,ranges=ranges),
                    **wlc_params)
    x_raw,y_raw = example.Separation,example.Force
    x0,model_x,model_y = FJC.fit_fjc_contour(x_raw,y_raw,**fit_dict)        
    print(x0)
    print(x0[0]/100)
    plt.plot(x_raw,y_raw)
    plt.plot(model_x,model_y)
    plt.show()
    
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
                                            cache_dir="./cache/",limit=3)
    ### XXX TODO
    # (1) correct for interference artifact
    # (2) get regions for WLC fit
    # (3) fit WLC to regions
    # (4) Invert WLC, determine dsDNA and ssDNA contour lengths at each force 
    region_fit_final = [8.1,8.9]
    region_fit_gc_rich = [7.6,7.8]
    region_fit_gc_poor = [7.37,7.49]
    # split the fecs...
    split_fecs = []
    for i,ex in enumerate(examples):
        split_fec = Analysis.zero_and_split_force_extension_curve(ex)
        retract = split_fec.retract
        split_fecs.append(split_fec)
    for s in split_fecs:
        retract = s.retract 
        max_fit_idx = np.argmax(retract.Force)
        plt.plot(retract.Time,retract.Force)
        plt.axvline(retract.Time[max_fit_idx])
        plt.show()
        continue
        for_wlc_fit = FEC_Util.slice_by_time(processed,*region_fit_final) 
        fit_polymer_model(for_wlc_fit)
    for i,ex in enumerate(examples):
        hairpin_plots(ex,filter_fraction=1e-3,out_path="./out/{:d}".format(i))
        
if __name__ == "__main__":
    run()
