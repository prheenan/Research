# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re,copy

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
   

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "./data_in/"
    data = IoUtilHao.read_and_cache_data_hao(in_dir,force=False)
    max_sep = 75e-9
    min_sep = 15e-9
    # get the slice we care about for each...
    for i,d in enumerate(data):
        sliced_fec = FEC_Util.slice_by_separation(d,min_sep,max_sep)
        fig = PlotUtilities.figure()
        plt.plot(d.Separation,d.Force,color='k',alpha=0.3)
        plt.plot(sliced_fec.Separation,sliced_fec.Force)
        plt.xlim([0,2*max_sep])
        plt.ylim([-50e-12,300e-12])
        PlotUtilities.savefig(fig,"out{:d}.png".format(i))

if __name__ == "__main__":
    run()
