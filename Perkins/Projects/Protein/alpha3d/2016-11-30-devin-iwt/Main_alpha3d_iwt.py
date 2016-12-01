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
    idx_end_of_unfolding = 7950
    unfold,refold = [],[]
    for d in RawData:
        # flip the sign, so force goes up 
        d.Force *= -1
        # get the unfolding and unfolds
        slice_unfolding = slice(0,idx_end_of_unfolding)
        unfold_tmp = FEC_Util.MakeTimeSepForceFromSlice(d,slice_unfolding)
        unfold_tmp.Separation -= np.percentile(unfold_tmp.Separation,5)
        unfold.append(unfold_tmp)
    # convert all the unfolding objects to IWT data
    IwtData = IWT_Util.ToIWTObjects(unfold)
    Landscape = InverseWeierstrass.FreeEnergyAtZeroForce(IwtData,NumBins=40)
    # get the titled landscape...
    Bounds = IWT_Util.BoundsObj(bounds_folded_nm=[0,10],
                                bounds_unfolded_nm=[0,10],
                                bounds_transition_nm=[0,10],
                                force_one_half_N=9e-12)
    tilt = IWT_Util.TiltedLandscape(Landscape,Bounds)
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(IwtData[-1].Extension*1e9,IwtData[-1].Force)
    plt.subplot(2,1,2)
    plt.plot(tilt.landscape_ext_nm,tilt.OffsetTilted)
    plt.show()
if __name__ == "__main__":
    run()
