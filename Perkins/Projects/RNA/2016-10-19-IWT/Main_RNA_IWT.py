# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")

from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import PlotUtilities

from scipy.signal import sawtooth

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
    for i in range(cycles):
        # XXX only get first half
        half = int((i+0.5)*cyclesize)
        fwd_slice = slice(i*cyclesize,half)
        rev_slice = slice(half,(i+1)*cyclesize)
        fwd = FEC_Util.MakeTimeSepForceFromSlice(Example,fwd_slice)
        rev = FEC_Util.MakeTimeSepForceFromSlice(Example,rev_slice)
        retracts.append(fwd)
        reverse.append(rev)
    half_cycle_size= fwd.Force.size
    up = np.linspace(z0,zf,half_cycle_size)
    down = np.linspace(zf,z0,half_cycle_size)
    updown = np.concatenate((up,down))
    cat_cyc = np.concatenate([updown for _ in range(cycles)])
    all_cycles = np.zeros(size)
    MaxSize = min(cat_cyc.size,all_cycles.size)
    all_cycles[:MaxSize] = cat_cyc[:MaxSize]
    # XXX just extrapolate end..
    all_cycles[MaxSize:] = cat_cyc[-1]
    IwtObjects = [InverseWeierstrass.\
        FEC_Pulling_Object(Time=o.Time,
                           Extension=o.Separation,
                           Force=o.Force,
                           SpringConstant=o.Meta.__dict__["K"],
                           Velocity=o.Meta.__dict__["RetractVelocity"],
                           ZFunc=lambda : up)
                  for o in retracts]
    InverseWeierstrass.SetAllWorkOfObjects(IwtObjects)
    # get the IWT
    Bins = [25,50,100,200,250,300,500,1000]
    for b in Bins:
        LandscapeObj =  InverseWeierstrass.\
                        FreeEnergyAtZeroForce(IwtObjects,
                                              NumBins=b)
        fig = PlotUtilities.figure(figsize=(8,8))
        IWT_Util.ForceExtensionHistograms(Example.Separation,
                                          Example.Force,
                                          nBins=b)
        PlotUtilities.savefig(fig,OutBase + "0_{:d}hist.pdf".format(b))
        fig = PlotUtilities.figure(figsize=(12,12))
        IWT_Util.EnergyLandscapePlot(LandscapeObj,FOneHalf=16e-12,
                                     ZoomBoundsMeters=[56e-9,68e-9],
                                     stiffness_pN_per_nm=80,
                                     NumPointsAround=int(b/10))
        PlotUtilities.savefig(fig,OutBase + "1_{:d}IWT.pdf".format(b))
    
if __name__ == "__main__":
    run()
