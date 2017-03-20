# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass

    
def hao_grouping_function(str_v):
    pattern = r"""
               (\D+) # type (eg: ext,force)
               (\d+)      # id (e.g. 1131)
               (\D+)        # anything else, who cares
               """
    match = re.match(pattern,str_v,re.VERBOSE)
    assert match is not None , "Whoops! Got a bad string: {:s}".format(str_v)
    ending,id,preamble = match.groups()
    # convert ext to sep
    ending = ending if ending != "ext" else "sep"
    return preamble,id,ending
    
def get_downsampled_data(downsample_n,force,cache_directory,
                         input_directory,limit):
    dict_kwargs = dict(cache_directory=cache_directory,
                       in_directory=input_directory,
                       grouping_function = hao_grouping_function,
                       limit=limit,force=force)
    data = FEC_Util.cache_ibw_directory(**dict_kwargs)  
    if (downsample_n > 1):
        get_down_slice = lambda x: slice(0,d.Force.size,downsample_n) 
        data = [FEC_Util.MakeTimeSepForceFromSlice(d,get_down_slice(d)) 
                for d in data]
    return data

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
    absolute_data_dir = base + \
        "4Patrick/Scratch/Tmp_Data_Scratch/hao-data-cache/"  
    downsample_n = 100
    fraction_for_vel = 0.1    
    limit = 10
    force_sample = False
    force= False
    cache_directory = "./"
    out_base = cache_directory    
    retracts = CheckpointUtilities.getCheckpoint("downsample.pkl",
                                                 get_downsampled_data,force,
                                                 downsample_n,force_sample,
                                                 cache_directory,
                                                 absolute_data_dir,limit=limit)
    # get just the retract regions
    # get after the surface invols
    retracts = [FEC_Util.GetFECPullingRegion(d,FlipSign=False,Correct=True) 
                for d in retracts]   
    # get the last time the data is above the median for each
    last_idx = []
    for r in retracts:
        force = r.Force
        x = r.Separation
        med = np.median(force)
        std = np.std(force)
        where_above = np.where(force > med + 3*std)[0]
        last_idx_above = where_above[-1]
        # determine where we go back the median after being way above items
        idx_arr = np.arange(force.size)
        where_back_below = np.where( (force <= med) & \
                                      (idx_arr > last_idx_above) )[0]
        last_idx_of_interest = where_back_below[0]                                 
        last_idx.append(last_idx_of_interest)
    # want to have a (small) amount after the last index, so use fudge
    fudge = lambda x: int(np.ceil(x.Force.size* 0.02))
    # slice the retract object to just where we care about
    retracts = [FEC_Util.MakeTimeSepForceFromSlice(d,slice(0,idx+fudge(d),1)) 
                for idx,d in zip(last_idx,retracts)]
    # comes in pN/nm, convert to N/m
    for r in retracts:
        r.Force *= 1e-12
        r.Separation *= 1e-9
    data_iwt = [IWT_Util.ToIWTObject(d) for d in retracts]
    # set all the effective velocities
    for d in data_iwt:
        IWT_Util.set_separation_velocity_by_first_frac(d,fraction_for_vel)
    # POST: they are all set. get the IWT 
    num_bins = 250
    LandscapeObj =  InverseWeierstrass.\
        FreeEnergyAtZeroForce(data_iwt,NumBins=num_bins)
    fig = PlotUtilities.figure(figsize=(12,12))
    IWT_Plot.plot_single_landscape(LandscapeObj,force_one_half_N=15e-12,
                                   add_meta_half=False,add_meta_free=False)  
    PlotUtilities.savefig(fig,out_base + "IWT.pdf")


if __name__ == "__main__":
    run()
