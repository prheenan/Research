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
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao

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
                                                 IoUtilHao.get_downsampled_data,
                                                 force,
                                                 downsample_n,force_sample,
                                                 cache_directory,
                                                 absolute_data_dir,limit=limit)
    data_iwt = [IoUtilHao.get_retract_pulling_region(d) for d in retracts]
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
