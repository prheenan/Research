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
from FitUtil.FitUtils.Python import FitUtil as pFitUtil
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

def ReadInAllFiles(FileNames,Limit):
    """
    Given a list of pxp files, reads them all into a list as 
    TimeSepForce Objcts

    Args:
        FileNames: List of .pxp full paths to data
        Limit: maximum number of curves to return
    """
    toRet = []
    for f in FileNames:
        toRet.extend(FEC_Util.ReadInData(f))
        # see if we are done
    # only return the limited number we want
    return toRet[:Limit]

def GetIWTObj(Base,FullNames,Force,Limit=150,
              PastZeroExt=60e-9,FilterToMeters=0.25e-9):
    """
    Given files, returns a list of IWT Objects for landscape reconstruction

    Args:
       Base: where to save the cache
       FullNames: ReadInAllFiles 
       Force: if true, force re-creation of the cache
       Limit: maximum number of curves to use 
       PastZeroExt: How much past the surface touchoff to analyze
       FilterToMeters: for determining the surface location, position resolution
       for a savitsky golay filter.
    Returns:
       list of IWT objects, list of TimeSepForce for the retraction, list of 
       TimeSepForce for just the desired amount past touchoff.
    """ 
    mObjs = pCheckUtil.getCheckpoint(Base + "cache.pkl",ReadInAllFiles,
                                     Force,FullNames,Limit)
    ApproachList,RetractList = FEC_Util.BreakUpIntoApproachAndRetract(mObjs)
    # filter all the retractions to the resolution specified,
    # based on the average velocity and the data sampling frequency.
    GetFilter = lambda x: max(3,
                              int((FilterToMeters/x.Velocity)*x.Frequency))
    Touchoff = [FEC_Util.GetFECPullingRegion(r,FilterPoints=GetFilter(r),
                                             MetersAfterTouchoff=PastZeroExt,
                                             Correct=True)
                for r in RetractList]
    # get the IWT transform objects
    IwtObjects = ToIWTObjects(Touchoff)
    return IwtObjects,RetractList,Touchoff

def GetObjectsAndIWT(Base,FullName,Force,NumBins=125,Limit=150):
    """
    Get the IWT and landscape associated with all the force extension curves

    Args:
        See GetIWTObj
    Returns:
        See GetIWTObj, plus it returns a LandscapeObj
    """
    IwtObjects,RetractList,Touchoff = GetIWTObj(Base,FullName,Force,
                                                Limit=Limit)
    # get the IWT
    LandscapeObj = InverseWeierstrass.FreeEnergyAtZeroForce(IwtObjects,
                                                            NumBins=NumBins)
    return IwtObjects,RetractList,Touchoff,LandscapeObj

def GetAllExtensionsAndForce(RetractList,Touchoff,IwtObjects,Base):
    toNano = lambda x : x * 1e9
    toPn = lambda x: x * 1e12
    NPlots = 3
    ext = []
    force = []
    MaxForce = max([np.max(t.Force) for t in Touchoff])
    MinForce = min([np.min(t.Force) for t in Touchoff])
    MaxX = max([t.Zsnsr[-1] for t in Touchoff])
    MaxWork = max([np.max(t.Work) for t in IwtObjects])
    # convert all the maxes to plot-friendly units
    MaxX_nm = MaxX * 1e9
    MaxWork_kbT = MaxWork/(4.1e-21)
    MaxForce_pN = MaxForce * 1e12
    MinForce_pN = MinForce * 1e12
    # get the limits
    ForceLim_pN = [MinForce_pN,MaxForce_pN]
    XLim_nm = [0,MaxX_nm]
    WorkLim_kbT = [0, MaxWork_kbT]
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
        plt.xlim(XLim_nm)
        # force bounds in pN
        plt.ylim(ForceLim_pN)
        pPlotUtil.lazyLabel("","Force (pN)","")
        plt.subplot(NPlots,1,3)
        plt.plot(toNano(Z),IwtObjects[i].Work/(4.1e-21),
                 alpha=0.3)
        plt.xlim(XLim_nm)
        plt.ylim(WorkLim_kbT)
        pPlotUtil.lazyLabel("Z stage Position (nm), relative to touchoff",
                            "Work (kbT)","")
        pPlotUtil.savefig(fig,Base + "{:d}.png".format(i))
        ext.extend(toNano(Z))
        force.extend(toPn(Touch.Force))
    return ext,force

def run():
    BaseDir = "/Volumes/group/4Patrick/Reports/" + \
           "2016_6_5_NearEquilibrium_Biotin_DOPC/IWT/"
    InBase = BaseDir + "In/"
    OutBase = BaseDir + "Out/"
    FullNames = [
    InBase +"2016-6-3-micah-1-part-per-million-biolevel-long-strept-coated.pxp",
    InBase +"2016-6-4-micah-1ppm-biolever-long-strept-saved-data.pxp",
    InBase +"2016-6-5-micah-1ppm-biolever-long-strept-saved-data.pxp"
    ]
    Limit = 100
    Stiffness_pN_per_nm = 4 # stiffness of the cantilever, pN/nm
    ForceReRead = False
    ForceRePlot = False
    IwtObjects,RetractList,Touchoff,LandscapeObj  = \
            pCheckUtil.getCheckpoint(OutBase + "IWT.pkl",GetObjectsAndIWT,
                                     ForceReRead,InBase,FullNames,ForceReRead,
                                     Limit=Limit)
    ext,force  = pCheckUtil.getCheckpoint(OutBase + "ExtAndForce.pkl",
                                          GetAllExtensionsAndForce,
                                          ForceRePlot,RetractList,Touchoff,
                                          IwtObjects,OutBase)
    # make a list of histograms for the plot
    fig = pPlotUtil.figure()
    ax = fig.add_subplot(111, projection='3d')
    H,xedges,yedges = np.histogram2d(force,ext,bins=20)
    NumX = xedges.size-1
    colors = ['r','g','b','y','m','k']
    NumColors = len(colors)
    for i in range(NumX):
        z = H[i,:]
        x = np.ones(z.size) * xedges[i]
        ax.bar(left=x,height=z,
               zs=yedges[:-1],zdir='y',alpha=0.7,color=colors[i%NumColors])
    pPlotUtil.lazyLabel("Force (pN)","Extension (nm)","",zlab="Count")
    pPlotUtil.savefig(fig,OutBase + "Hist2d.png")
    toNano = lambda x : x * 1e9
    toPn = lambda x: x * 1e12
    nBins = 50
    # make a heat map, essentially
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    counts, xedges, yedges, Image = plt.hist2d(ext, force,
                                               bins=nBins,cmap='afmhot')
    x_bins = xedges[:-1]
    y_bins = yedges[:-1]
    bindiff = np.median(np.diff(x_bins))
    for i in range(x_bins.size):
        N = sum(counts[i,:])
        average = (sum(counts[i,:] * y_bins))/N
        label = "Avg binned force" if i == 0 else ""
        plt.plot(x_bins[i]+bindiff/2,average,'go',label=label)
    pPlotUtil.lazyLabel("Tip-Bilayer Separation [nm]",
                        "Force [pN]",
                        "Two-Dimensional Force-Separation Histogram",
                        frameon=True)
    cbar = plt.colorbar()
    cbar.set_label('# in (Force,Separation) Bin', labelpad=10,rotation=270)
    pPlotUtil.savefig(fig,OutBase + "HeatMap.png")
    fig = pPlotUtil.figure(figsize=(8,12))
    plt.subplot(3,1,1)
    NanoExt =toNano(LandscapeObj.Extensions)
    FreeEnergyEq = LandscapeObj.EnergyLandscape
    plt.plot(NanoExt,FreeEnergyEq * LandscapeObj.Beta)
    pPlotUtil.lazyLabel("","G0",
                        "Reconstructed Energy Landscape for lipid Bilayer")
    plt.subplot(3,1,2)
    FOneHalf = 8e-12
    TiltedEnergy = (FreeEnergyEq - LandscapeObj.Extensions * FOneHalf)
    TiltedEnergy -= (TiltedEnergy[0])
    EnergyByBeta = TiltedEnergy * LandscapeObj.Beta
    plt.plot(NanoExt,EnergyByBeta)
    pPlotUtil.lazyLabel("Molecular Extension (nm)","G at F-1/2 (kT)",
                        "")
    plt.subplot(3,1,3)
    # zoom on in a specific reigon, just eyeballing
    # the bounds based on the 'full' landscape plot
    WhereZoom = np.where( (NanoExt > 22) &
                          (NanoExt < 30))
    EnergyZoom = EnergyByBeta[WhereZoom]
    ExtZoom = NanoExt[WhereZoom]
    # Zero out everything
    MinIdx = np.argmin(EnergyZoom)
    EnergyZoom -= np.min(EnergyZoom)
    ExtZoom -= ExtZoom[MinIdx]
    # fit a parabola to the bottom of the well
    NumPointsAround =  4
    IdxSlice = slice(MinIdx-NumPointsAround,MinIdx+NumPointsAround)
    FitX = ExtZoom[IdxSlice]
    FitY = EnergyZoom[IdxSlice]
    coeffs = np.polyfit(FitX,FitY,deg=2)
    # remove the linear term, which is the second
    coeffs[1] = 0
    xinterp = np.linspace(FitX[0],FitX[-1])
    vals = np.polyval(coeffs,x=xinterp)
    curvature_kbT_per_nm = coeffs[0]
    # note: the well is in kT/nm, so we convert to pN/nm
    # by multuplying by 4.1
    curvature_pN_per_nm = curvature_kbT_per_nm * 4.1
    plt.plot(ExtZoom,EnergyZoom)
    plt.plot(xinterp,vals,color='g',linewidth=4.0,linestyle='--',
             label="Landscape parabola: ({:.1f} pN/nm)".\
             format(curvature_pN_per_nm))
    plt.plot(xinterp,Stiffness_pN_per_nm*xinterp**2+max(vals),
             label="Cantilever Response ({:.1f} pN/nm)".\
             format(Stiffness_pN_per_nm))
    plt.ylim([0,max(EnergyZoom)])
    pPlotUtil.lazyLabel("Molecular Extension (nm)","G at F-1/2 (kT)",
                        "",frameon=True)
    pPlotUtil.savefig(fig,OutBase + "IWT.png")

if __name__ == "__main__":
    run()
