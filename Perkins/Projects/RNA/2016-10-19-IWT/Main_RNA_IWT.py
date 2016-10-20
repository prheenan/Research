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
    z0_Separation = 59e-9
    zf_Separation = 75e-9
    z0_ZSnsr = 57e-9
    zf_ZSnsr = 70e-9
    # choose the bounds
    z0 = z0_Separation
    zf = zf_Separation
    cycles = 16
    size = Example.ZSnsr.size
    cyclesize = int(size/cycles)
    # XXX break up by cycle
    for i in range(cycles):
        m_slice = slice(i*cyclesize,(i+1)*cyclesize)
        plt.plot(Example.Force[m_slice])
        plt.show()
    half_cycle_size= cyclesize/2
    up = np.linspace(z0,zf,half_cycle_size)
    down = np.linspace(zf,z0,half_cycle_size)
    updown = np.concatenate((up,down))
    cat_cyc = np.concatenate([updown for _ in range(cycles)])
    all_cycles = np.zeros(size)
    MaxSize = min(cat_cyc.size,all_cycles.size)
    all_cycles[:MaxSize] = cat_cyc[:MaxSize]
    print(all_cycles.size-cat_cyc.size)
    # XXX just extrapolate end..
    all_cycles[MaxSize:] = cat_cyc[-1]
    plt.plot(Example.Separation)
    plt.plot(RawData[0].Separation,color='k',alpha=0.3)
    plt.plot(cat_cyc)
    plt.show()
    IwtObjects = [InverseWeierstrass.\
        FEC_Pulling_Object(Time=o.Time,
                           Extension=o.Separation,
                           Force=o.Force,
                           SpringConstant=o.Meta.__dict__["K"],
                           Velocity=o.Meta.__dict__["RetractVelocity"],
                           ZFunc=lambda : all_cycles)
                  for o in FilteredData]
    InverseWeierstrass.SetAllWorkOfObjects(IwtObjects)
    # get the IWT
    Bins = [200,250,300]
    for b in Bins:
        LandscapeObj =  InverseWeierstrass.\
                        FreeEnergyAtZeroForce(IwtObjects,
                                              NumBins=b)
        fig = PlotUtilities.figure(figsize=(8,8))
        IWT_Util.ForceExtensionHistograms(Example.Separation,
                                          Example.Force,
                                          nBins=b)
        PlotUtilities.savefig(fig,OutBase + "0_{:d}hist.pdf".format(b))
        fig = PlotUtilities.figure(figsize=(8,8))
        IWT_Util.EnergyLandscapePlot(LandscapeObj,FOneHalf=16e-12)
        PlotUtilities.savefig(fig,OutBase + "1_{:d}IWT.pdf".format(b))
    
if __name__ == "__main__":
    run()
