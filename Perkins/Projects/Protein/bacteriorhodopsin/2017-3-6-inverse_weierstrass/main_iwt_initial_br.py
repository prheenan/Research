# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao
from FitUtil.WormLikeChain.Python.Code import WLC 
import cProfile 
   
def get_first_peak_slices(absolute_data_dir):
    downsample_n = 1
    fraction_for_vel = 0.1    
    limit = 5
    force_sample = False
    force= False
    cache_directory = "./"
    out_base = cache_directory    
    retracts = IoUtilHao.get_downsampled_data(downsample_n,force_sample,
                                              cache_directory,
                                              absolute_data_dir,limit=limit)
    data_fec = [IoUtilHao.get_retract_pulling_region(d) for d in retracts]
    # get just the first X nm
    max_sep = 30e-9
    slice_func = lambda x: slice(0,np.where(x.Separation >= max_sep)[0][0],1)
    just_first = [FEC_Util.MakeTimeSepForceFromSlice(r,slice_func(r)) 
                  for r in data_fec]
    # get up until the maximum in this region
    slice_max = lambda x: slice(0,np.argmax(x.Force),1)
    until_max = [FEC_Util.MakeTimeSepForceFromSlice(r,slice_max(r)) 
                 for r in just_first]
    data_iwt = []
    for r in until_max:
        tmp_iwt = IWT_Util.ToIWTObject(r)
        # set all the effective velocities
        IWT_Util.set_separation_velocity_by_first_frac(tmp_iwt,fraction_for_vel)
        data_iwt.append(tmp_iwt)
    return data_iwt
    
def get_contour_lengths(data_iwt):
    contour_lengths = []
    for reference in data_iwt:
        max_ext = max(reference.Extension)
        brute_dict = dict(ranges=[ [max_ext/2,max_ext*2] ],Ns=15)
        kwargs_fit = dict(kbT = 4.11e-21,
                          Lp=0.3e-9,
                          K0=1000e-12)
        sep,force = reference.Separation,reference.Force                      
        x0,y = WLC.fit(sep,force,brute_dict=brute_dict,
                       **kwargs_fit)
        contour_lengths.append(x0)
    return np.concatenate(contour_lengths)
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = FEC_Util.default_data_root()
    # XXX use Haos...
    absolute_data_dir = base + "4Patrick/BRFEC/FEC_3ums/"
    out_base = "./"
    # get the region near the first peak
    data_iwt = CheckpointUtilities.getCheckpoint("./peaks.pkl",
                                                 get_first_peak_slices,False,
                                                 absolute_data_dir)
    # get the contour lengths                                                 
    contour_lengths = \
        CheckpointUtilities.getCheckpoint("./contour.pkl",
                                          get_contour_lengths,False,data_iwt)
    
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(data_iwt)
    PlotUtilities.savefig(fig,"heatmap.png")                 
    # POST: they are all set. get the IWT 
    num_bins = 250
    pr = cProfile.Profile()
    pr.enable()
    LandscapeObj =  InverseWeierstrass.\
        FreeEnergyAtZeroForce(data_iwt,NumBins=num_bins)
    pr.disable()
    pr.print_stats(sort='time') 
    force_N = [1e-12,5e-12,10e-12,20e-12,40e-12,100e-12,250e-12,500e-12]
    for i,f in enumerate(force_N):
        fig = PlotUtilities.figure(figsize=(12,12))

        IWT_Plot.plot_single_landscape(LandscapeObj,f_one_half_N=f,
                                       add_meta_half=False,add_meta_free=False) 
        out_name= out_base + "IWT{:d}_{:.1g}.png".format(i,f*1e12)   
        PlotUtilities.savefig(fig,out_name)


if __name__ == "__main__":
    run()
