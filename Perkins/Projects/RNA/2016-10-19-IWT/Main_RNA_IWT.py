# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")

from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import PlotUtilities

from scipy.signal import sawtooth
# file-wide bounds on our state extensions in nm
bounds_folded_nm = [56,61]
bounds_unfolded_nm = [64,69]
bounds_transition_nm = [61.75,63.5]
FOneHalf_N = 14.2e-12
    
def run():
    """
    """
    Base = "./"
    OutBase = Base + "out/"
    InFiles = [Base + "RNA_Ramp_20nmps.pxp"]
    RawData = IWT_Util.ReadInAllFiles(InFiles,Limit=50)
    FilterPoints = 100
    FilteredData = [FEC_Util.GetFilteredForce(r,FilterPoints)
                    for r in RawData]
    Example = FilteredData[0]
    z0_Separation = 60e-9
    zf_Separation = 73e-9
    z0_ZSnsr = 57e-9
    zf_ZSnsr = 70e-9
    # choose the bounds
    z0 = z0_Separation
    zf = zf_Separation
    cycles = 16
    size = Example.ZSnsr.size
    cyclesize = int(size/cycles)
    retracts,reverse = [],[]
    bounds = IWT_Plot.BoundsObj(bounds_folded_nm,bounds_unfolded_nm,
                                bounds_transition_nm,FOneHalf_N)
    for i in range(cycles):
        # XXX only get first half
        half = int((i+0.5)*cyclesize)
        fwd_slice = slice(i*cyclesize,half)
        rev_slice = slice(half,(i+1)*cyclesize)
        fwd = FEC_Util.MakeTimeSepForceFromSlice(Example,fwd_slice)
        rev = FEC_Util.MakeTimeSepForceFromSlice(Example,rev_slice)
        retracts.append(fwd)
        reverse.append(rev)
    fwd_cycle_size= fwd.Force.size
    rev_cycle_size = rev.Force.size
    up = np.linspace(z0,zf,fwd_cycle_size)
    down = np.linspace(zf,z0,rev_cycle_size)
    updown = np.concatenate((up,down))
    cat_cyc = np.concatenate([updown for _ in range(cycles)])
    all_cycles = np.zeros(size)
    MaxSize = min(cat_cyc.size,all_cycles.size)
    all_cycles[:MaxSize] = cat_cyc[:MaxSize]
    # XXX just extrapolate end..
    all_cycles[MaxSize:] = cat_cyc[-1]
    UnfoldObj = [IWT_Util.RobTimeSepForceToIWT(o,ZFunc=(lambda: up))
                 for o in retracts]
    RefoldObj = [IWT_Util.RobTimeSepForceToIWT(o,ZFunc=(lambda: down))
                 for o in reverse]
    LandscapeObj =  InverseWeierstrass.\
        FreeEnergyAtZeroForce(UnfoldObj,NumBins=40,RefoldingObjs=RefoldObj)
    for idx in range(min(len(UnfoldObj),RefoldObj)):
         IWT_Plot.TomPlot(LandscapeObj,OutBase+ str(idx) + "_",UnfoldObj,
                          RefoldObj,idx=idx,bounds=bounds)
    IWT_Plot.InTheWeedsPlot(OutBase ,UnfoldObj,RefoldObj,Example,
                            Bins = [50,75,100,150],bounds=bounds)


if __name__ == "__main__":
    run()
