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

import cProfile 
   


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
    downsample_n = 1
    fraction_for_vel = 0.1    
    limit = 30
    force_sample = False
    force= False
    cache_directory = "./"
    out_base = cache_directory    
    retracts = IoUtilHao.get_downsampled_data(downsample_n,force_sample,
                                              cache_directory,
                                              absolute_data_dir,limit=limit)
    data_iwt = [IoUtilHao.get_retract_pulling_region(d) for d in retracts]
    # POST: they are all set. get the IWT 
    num_bins = 250
    pr = cProfile.Profile()
    pr.enable()
    LandscapeObj =  InverseWeierstrass.\
        FreeEnergyAtZeroForce(data_iwt,NumBins=num_bins)
    pr.disable()
    pr.print_stats(sort='time')
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(data_iwt)
    PlotUtilities.savefig(fig,"heatmap.png")    
    force_N = [1e-12,5e-12,10e-12,20e-12,40e-12,100e-12,250e-12,500e-12]
    for i,f in enumerate(force_N):
        fig = PlotUtilities.figure(figsize=(12,12))

        IWT_Plot.plot_single_landscape(LandscapeObj,f_one_half_N=f,
                                       add_meta_half=False,add_meta_free=False) 
        out_name= out_base + "IWT{:d}_{:.1g}.png".format(i,f*1e12)   
        PlotUtilities.savefig(fig,out_name)


if __name__ == "__main__":
    run()
