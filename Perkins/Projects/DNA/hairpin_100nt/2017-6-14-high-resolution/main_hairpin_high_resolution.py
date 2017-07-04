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
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
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
    """
    Fits the polymer model de jour to example
    
    Args:
        example: the timesepforce to fit to
    Returns:
        x0, the parameters of the fit 
    """
    wlc_params = dict(K0=2000e-12,kbT=4.1e-21)
    ranges = [(10e-9,90e-9),(0.1e-9,1e-9)]
    fit_dict = dict(brute_dict=dict(Ns=30,ranges=ranges),
                    **wlc_params)
    x_raw,y_raw = example.Separation,example.Force
    x0,model_x,model_y = FJC.fit_fjc_contour(x_raw,y_raw,**fit_dict)        
    return x0
    
def get_polymer_coefficients(split_fecs,working_distance_nm):
    """
    returns the polymer coefficients for the given examples, as well as the 
    'valid' indices (ie: if we couldn't fit fot some reason, that idx wont be 
    there)
    
    Args;
        examples:
        working_distance_nm: currently, we assume     
    Returns:
        tuple of <coefficients, valid indices>
    """
    coeffs,idx_valid = [],[]
    for i,s in enumerate(split_fecs):
        retract = s.retract 
        max_fit_idx = np.argmax(retract.Force)
        # find the first idx we are <= the working distance
        target_sep = (retract.Separation[max_fit_idx] - \
                      working_distance_nm * 1e-9)
        idx = np.where(retract.Separation <= target_sep)[0]
        if (idx.size == 0):
            continue
        # POST: at least one thing to fit 
        fit_slice = slice(idx[-1],max_fit_idx,1)
        retract_fit = FEC_Util.MakeTimeSepForceFromSlice(retract,fit_slice)
        plt.plot(retract_fit.Separation,retract_fit.Force)
        x0 = fit_polymer_model(retract_fit)
        # offset the retract by the contour length (first paramter)
        coeffs.append(x0)
        idx_valid.append(i)
    return coeffs,idx_valid
             
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    abs_dir = "./"
    cache_dir = "./cache/"
    examples = FEC_Util.\
        cache_individual_waves_in_directory(pxp_dir=abs_dir,force=False,
                                            cache_dir=cache_dir,limit=20)
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
    # align them by the contour lengths (assuming we have at least xnm to 
    # work with 
    working_distance_nm = 30
    coeffs,idx = CheckpointUtilities.\
        getCheckpoint("./polymer.pkl",get_polymer_coefficients,False,
                      split_fecs,working_distance_nm)
    # get the split fecs we could actually fit
    good_splits = [split_fecs[i] for i in idx]
    # align all the retracts by the contour lenghts
    contour_L0 = [c[0] for c in coeffs]
    arbitrary_offset = 75e-9
    for L0,split_fec in zip(contour_L0,good_splits):
        split_fec.retract.Separation -= L0
        split_fec.retract.Separation += arbitrary_offset
    # POST: retracts are all aligned. 
    # make a heat map 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec([r.retract for r in good_splits])
    PlotUtilities.title("FEC Heat map, aligned by L0, N={:d}".format(len(good_splits)))
    PlotUtilities.savefig(fig,"./heat.png")
    for i,ex in enumerate(examples):
        hairpin_plots(ex,filter_fraction=1e-3,out_path="./out/{:d}".format(i))
    # XXX plot the data with the fit of the WLC aligned ontop 
        
if __name__ == "__main__":
    run()
