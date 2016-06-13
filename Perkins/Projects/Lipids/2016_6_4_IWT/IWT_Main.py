# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D



sys.path.append("../../../../..")
# XXX move to utility?

from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
from GeneralUtil.python import PlotUtilities as pPlotUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot

def ToIWTObjects(TimeSepForceObjects):
    """
    Converts TimeSepForceObjects to InverseWeierstrass objects

    Args:
        TimeSepForceObjects: list of TimeSepForceObjects to transform
    """
    Objs = [InverseWeierstrass.\
            FEC_Pulling_Object(Time=o.Time,
                               Extension=o.Zsnsr,
                               Force=o.Force,
                               SpringConstant=o.SpringConstant,
                               Velocity=o.Velocity)
            for o in TimeSepForceObjects]
    InverseWeierstrass.SetAllWorkOfObjects(Objs)
    return Objs

def GetIWTObj(Base,FullName):
    Limit=150
    mObjs = pCheckUtil.getCheckpoint(Base + "cache.pkl",FEC_Util.ReadInData,
                                     False,FullName,Limit)
    ApproachList,RetractList = FEC_Util.BreakUpIntoApproachAndRetract(mObjs)
    # filter all the retractions
    PastZeroExt = 8e-9
    FilterPoints = 30
    # get just after the touchoff
    FilterToMeters = 0.25e-9
    GetFilter = lambda x: max(3,
                              int((FilterToMeters/x.Velocity)*x.Frequency))
    Touchoff = [FEC_Util.GetFECPullingRegion(r,FilterPoints=GetFilter(r),
                                             MetersAfterTouchoff=PastZeroExt)
                for r in RetractList]
    # get the IWT transform objects
    IwtObjects = ToIWTObjects(Touchoff)
    return IwtObjects,RetractList,Touchoff

def GetObjectsAndIWT(Base,FullName):
    IwtObjects,RetractList,Touchoff = GetIWTObj(Base,FullName)
    # get the IWT
    LandscapeObj = InverseWeierstrass.FreeEnergyAtZeroForce(IwtObjects,
                                                            NumBins=100)
    return IwtObjects,RetractList,Touchoff,LandscapeObj

def GetAllExtensionsAndForce(RetractList,Touchoff,IwtObjects,Base):
    toNano = lambda x : x * 1e9
    toPn = lambda x: x * 1e12
    NPlots = 3
    xlim_nm = [0,10]
    ext = []
    force = []
    for i,(Retract,Touch) in enumerate(zip(RetractList,Touchoff)):
        fig = pPlotUtil.figure(figsize=(8,12))
        ForZ = Retract
        RetractZ = Retract.Zsnsr
        RetractZ -= np.min(RetractZ)
        plt.subplot(NPlots,1,1)
        # normalize and flip the force, XXX move to utility...
        ForceRetractPlot = toPn(Retract.Force)
        ForceRetractPlot *= -1
        N = ForceRetractPlot.size
        fraction = 0.2
        ForceRetractPlot -= np.median(ForceRetractPlot[-int(fraction*N):])
        plt.plot(toNano(RetractZ),ForceRetractPlot,alpha=0.3)
        pPlotUtil.lazyLabel("Z stage Position (nm), Absolute","Force (pN)",
                            "Determining Work and Force for FEC")
        plt.subplot(NPlots,1,2)
        Z = Touch.Zsnsr
        plt.plot(toNano(Z),toPn(Touch.Force),alpha=0.3)
        plt.xlim(xlim_nm)
        # force bounds in pN
        plt.ylim([-25,30])
        pPlotUtil.lazyLabel("","Force (pN)","")
        plt.subplot(NPlots,1,3)
        plt.plot(toNano(Z),IwtObjects[i].Work/(4.1e-21),
                 alpha=0.3)
        plt.xlim(xlim_nm)
        plt.ylim([-2,20])
        pPlotUtil.lazyLabel("Z stage Position (nm), relative to touchoff",
                            "Work (kbT)","")
        pPlotUtil.savefig(fig,Base + "{:d}.png".format(i))
        ext.extend(toNano(Z))
        force.extend(toPn(Touch.Force))
    return ext,force

def run():
    Base = "/Users/patrickheenan/Documents/education/boulder_files/" +\
           "rotations_year_1/3_perkins/reports/2016_Bio-DOPC-Energy-Landscape/"
    FullName = Base + "prh_cleaned_2016-6-4-micah-1ppm-biolever-long-" +\
               "strept-saved-data.pxp"
    IwtObjects,RetractList,Touchoff,LandscapeObj  = \
            pCheckUtil.getCheckpoint(Base + "IWT.pkl",GetObjectsAndIWT,
                                     False,Base,FullName)
    ext,force  = pCheckUtil.getCheckpoint(Base + "ExtAndForce.pkl",
                                          GetAllExtensionsAndForce,
                                          True,RetractList,Touchoff,IwtObjects,
                                          Base)
    # bin the force by extensions
    fig = pPlotUtil.figure()
    ax = fig.add_subplot(111, projection='3d')
    H,xedges,yedges = np.histogram2d(force,ext,bins=20)
    NumX = xedges.size-1
    colors = ['r','g','b','y','m']
    NumColors = len(colors)
    for i in range(NumX):
        z = H[i,:]
        x = np.ones(z.size) * xedges[i]
        ax.bar(left=x,height=z,
               zs=yedges[:-1],zdir='y',alpha=0.7,color=colors[i%NumColors])
    pPlotUtil.lazyLabel("Force (pN)","Extension (nm)","",zlab="Count")
    pPlotUtil.savefig(fig,Base + "Hist2d.png")
    toNano = lambda x : x * 1e9
    toPn = lambda x: x * 1e12
    nBins = 40
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    cax = plt.hist2d(ext, force, bins=nBins,cmap='afmhot')
    pPlotUtil.lazyLabel("Tip-Bilayer Separation [nm]",
                        "Force [pN]",
                        "Two-Dimensional Force-Separation Histogram")
    cbar = plt.colorbar()
    cbar.set_label('# in (Force,Separation) Bin', labelpad=10,rotation=270)
    pPlotUtil.savefig(fig,Base + "HeatMap.png")
    fig = pPlotUtil.figure(figsize=(8,12))
    xlim = [-1,7]
    plt.subplot(2,1,1)
    NanoExt =toNano(LandscapeObj.Extensions)
    FreeEnergyEq = LandscapeObj.EnergyLandscape
    plt.plot(NanoExt,FreeEnergyEq * LandscapeObj.Beta)
    plt.xlim(xlim)
    pPlotUtil.lazyLabel("","G0",
                        "DeltaG for lipid rupture is approximately 1kT")
    plt.subplot(2,1,2)
    FOneHalf = 8e-12
    TiltedEnergy = (FreeEnergyEq - LandscapeObj.Extensions * FOneHalf)
    TiltedEnergy -= (TiltedEnergy[0])
    plt.plot(NanoExt,TiltedEnergy * LandscapeObj.Beta)
    plt.xlim(xlim)
    pPlotUtil.lazyLabel("Molecular Extension (nm)","G at F-1/2 (kT)",
                        "")
    fig.savefig(Base + "IWT.png")

if __name__ == "__main__":
    run()
