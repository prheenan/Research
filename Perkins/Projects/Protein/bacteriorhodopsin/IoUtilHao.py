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
from Research.Personal.EventDetection.Util import Analysis
from Research.Personal.EventDetection._2SplineEventDetector import _no_event

    
def hao_grouping_function(str_v):
    pattern = r"""
               (ext|force)     # type (eg: ext,force)
               \D+?            # non-greedy non digits (e.g. "extsmooth")
               (\d+)           # id (e.g. 1131)
               (\D+)           # anything else, who cares
               """
    match = re.match(pattern,str_v,re.VERBOSE)
    assert match is not None , "Whoops! Got a bad string: {:s}".format(str_v)
    ending,id,preamble = match.groups()
    # convert ext to sep
    ending = ending if ending != "ext" else "sep"
    return preamble,id,ending
    
def read_and_cache_data_hao(input_directory,cache_directory="./cache_in/",
                            force=True,limit=None,**kwargs):
    """
    Reads hao-style ibw files from input_directory, caching them as TimeSepForce
    to cache_directory 
    
    Args:
        see cache_ibw_directory
    Returns;
        list of TimeSepForce objects 
    """    
    dict_kwargs = dict(cache_directory=cache_directory,
                       in_directory=input_directory,
                       grouping_function = hao_grouping_function,
                       limit=limit,force=force)
    data = FEC_Util.cache_ibw_directory(**dict_kwargs)  
    for d in data:
        d.Force *= 1e-12
        d.Separation *= 1e-9
    return data
    
def get_downsampled_data(downsample_n,force,cache_directory,
                         input_directory,limit):
    data = read_data_hao(**kwargs)
    if (downsample_n > 1):
        get_down_slice = lambda x: slice(0,d.Force.size,downsample_n) 
        data = [FEC_Util.MakeTimeSepForceFromSlice(d,get_down_slice(d)) 
                for d in data]
    return data
    
def get_retract_pulling_region(d,fraction_fudge=0.02,
                               n_std=3,fraction_zero_retract=0.05,zero=False):
    # get just the retract regions
    # get after the surface invols
    r = FEC_Util.GetFECPullingRegion(d,FlipSign=False,Correct=True)
    # get the last time the data is above the median for each
    force = r.Force
    x = r.Separation
    med = np.median(force)
    std = np.std(force)
    where_above = np.where(force > med + n_std*std)[0]
    last_idx_above = where_above[-1]
    # determine where we go back the median after being way above items
    idx_arr = np.arange(force.size)
    where_back_below = np.where( (force <= med) & \
                                  (idx_arr > last_idx_above) )[0]
    last_idx = where_back_below[0]                                 
    # want to have a (small) amount after the last index, so use fudge
    fudge = lambda x: int(np.ceil(x.Force.size* fraction_fudge))
    # slice the retract object to just where we care about
    idx_f = last_idx+fudge(r)
    retract = FEC_Util.MakeTimeSepForceFromSlice(r,slice(0,idx_f,1))
    effective_approach_feather = \
        FEC_Util.MakeTimeSepForceFromSlice(r,slice(idx_f,None,1))
    # comes in pN/nm, convert to N/m
    retract.Force *= 1e-12
    retract.Separation *= 1e-9
    # zero by the bottom P%...
    if (zero):
        retract.Separation -= np.percentile(retract.Separation,
                                            fraction_zero_retract)
    return retract


if __name__ == "__main__":
    run()
