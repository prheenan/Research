# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from IgorUtil.PythonAdapter import PxpLoader
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    Base = "./"
    OutBase = Base + "out/"
    InFiles = [Base + "ForPatrick.pxp"]
    RawData = IWT_Util.\
              ReadInAllFiles(InFiles,Limit=50,
                             ValidFunc=PxpLoader.valid_fec_allow_endings)
    # get the start/ends of the re-folding and unfolding portions
    # hard coded constant for now...
    # XXX for re-folding, need to add in schedule
    # XXX ensure properly zeroed?
    idx_end_of_unfolding = 8100
    unfold,refold = [],[]
    for d in RawData:
        # flip the sign, so force goes up 
        d.Force *= -1
        # get the unfolding and unfolds
        slice_unfolding = slice(0,idx_end_of_unfolding)
        unfold_tmp = FEC_Util.MakeTimeSepForceFromSlice(d,slice_unfolding)
        slice_folding = slice(idx_end_of_unfolding,idx_end_of_unfolding*2)
        fold_tmp = FEC_Util.MakeTimeSepForceFromSlice(d,slice_folding)
        # arbitrarily assign the minimum separaiton to the lower 5%.
        MinV = np.percentile(unfold_tmp.Separation,5)
        unfold_tmp.Separation -= MinV
        fold_tmp.Separation -= MinV
        unfold.append(unfold_tmp)
        refold.append(fold_tmp)
    # convert all the unfolding objects to IWT data
    IwtData = IWT_Util.ToIWTObjects(unfold)
    IwtData_fold = IWT_Util.ToIWTObjects(refold)
    # switch the velocities of all the folding objects..
    for o in IwtData_fold:
        o.Velocity *= -1
    # get the titled landscape...
    Bounds = IWT_Util.BoundsObj(bounds_folded_nm=[20,30],
                                bounds_transition_nm=[26,35],
                                bounds_unfolded_nm=[32,40],
                                force_one_half_N=13e-12)
    IWT_Plot.InTheWeedsPlot(OutBase="./out/",
                            UnfoldObj=IwtData,
                            bounds=Bounds,Bins=[40,60,80],
                            max_landscape_kT=None,
                            min_landscape_kT=None)

if __name__ == "__main__":
    run()
