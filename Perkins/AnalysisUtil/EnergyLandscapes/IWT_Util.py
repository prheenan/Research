# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt


from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
from GeneralUtil.python import PlotUtilities
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

def toNano(x):
    return x * 1e9


def toPn(x):
    return x * 1e12


def GetAllExtensionsAndForceAndPlot(RetractList,Touchoff,IwtObjects,Base):
    """
    Returns all the normalized forces and extnsions, after plotting them
    in the location given by base

    Args:
        RetractList: list of TimeSepForce objects with just the retract
        Touchoff: list of TimeSepForceObjects *zero-offset* in x and y
        IwtObjects: Touchoff, converted to IWT
        Base: where to save
    Returns:
        tuple of <stage Z,force>
    """
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
        fig = PlotUtilities.figure(figsize=(8,12))
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
        PlotUtilities.lazyLabel("Z stage Position (nm), Absolute","Force (pN)",
                            "Determining Work and Force for FEC")
        plt.subplot(NPlots,1,2)
        Z = Touch.Zsnsr
        plt.plot(toNano(Z),toPn(Touch.Force),alpha=0.3)
        plt.xlim(XLim_nm)
        # force bounds in pN
        plt.ylim(ForceLim_pN)
        PlotUtilities.lazyLabel("","Force (pN)","")
        plt.subplot(NPlots,1,3)
        plt.plot(toNano(Z),IwtObjects[i].Work/(4.1e-21),
                 alpha=0.3)
        plt.xlim(XLim_nm)
        plt.ylim(WorkLim_kbT)
        PlotUtilities.lazyLabel("Z stage Position (nm), relative to touchoff",
                            "Work (kbT)","")
        PlotUtilities.savefig(fig,Base + "{:d}.png".format(i))
        ext.extend(toNano(Z))
        force.extend(toPn(Touch.Force))
    return ext,force

def ForceExtensionHistograms(ext,force,nBins=100):
    """
    Makes a 2-d force histogram (ext,force)
    
    Args:
        ext: list of extensions
        force: list of forces
        nBins: how many bins to use
    """
    # make a heat map, essentially
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
    PlotUtilities.lazyLabel("Tip-Bilayer Separation [nm]",
                        "Force [pN]",
                        "Two-Dimensional Force-Separation Histogram",
                        frameon=True)
    cbar = plt.colorbar()
    cbar.set_label('# in (Force,Separation) Bin', labelpad=10,rotation=270)


def EnergyLandscapePlot(LandscapeObj,FOneHalf=8e-12,
                        ZoomBoundsMeters=[22e-9,30e-9],
                        NumPointsAround=4,
                        stiffness_pN_per_nm=4):
    """
    Plots the enegry landscape  and tilted landscape

    Args:
        LandscapeObj: return from InverseWeierstrass.FreeEnergyAtZeroForce(
        FOneHalf: what to tilt by
    """
    plt.subplot(3,1,1)
    NanoExt =toNano(LandscapeObj.Extensions)
    FreeEnergyEq = LandscapeObj.EnergyLandscape
    plt.plot(NanoExt,FreeEnergyEq * LandscapeObj.Beta)
    PlotUtilities.lazyLabel("","G0",
                        "Reconstructed Energy Landscape for lipid Bilayer")
    plt.subplot(3,1,2)
    TiltedEnergy = (FreeEnergyEq - LandscapeObj.Extensions * FOneHalf)
    TiltedEnergy -= (TiltedEnergy[0])
    EnergyByBeta = TiltedEnergy * LandscapeObj.Beta
    plt.plot(NanoExt,EnergyByBeta)
    PlotUtilities.lazyLabel("Molecular Extension (nm)","G at F-1/2 (kT)",
                            "")
    plt.subplot(3,1,3)
    # zoom on in a specific reigon, just eyeballing
    # the bounds based on the 'full' landscape plot
    ZoomNm = np.array(ZoomBoundsMeters)* 1e9
    WhereZoom = np.where( (NanoExt > min(ZoomNm)) &
                          (NanoExt < max(ZoomNm)))
    if (WhereZoom[0].size > 0):
        EnergyZoom = EnergyByBeta[WhereZoom]
        ExtZoom = NanoExt[WhereZoom]
        # Zero out everything
        MinIdx = np.argmin(EnergyZoom)
        EnergyZoom -= np.min(EnergyZoom)
        ExtZoom -= ExtZoom[MinIdx]
        # fit a parabola to the bottom of the well
        IdxStart = max(0,MinIdx-NumPointsAround)
        IdxEnd = min(ExtZoom.size,MinIdx+NumPointsAround)
        IdxSlice = slice(IdxStart,IdxEnd)
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
        plt.plot(xinterp,stiffness_pN_per_nm*xinterp**2+max(vals),
                 label="Cantilever Response ({:.1f} pN/nm)".\
                 format(stiffness_pN_per_nm))
        plt.ylim([0,max(EnergyZoom)])
        PlotUtilities.lazyLabel("Molecular Extension (nm)","G at F-1/2 (kT)",
                                "",frameon=True)




