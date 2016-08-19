# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis.DataCorrection.\
    CorrectionMethods import CorrectForcePullByMetaInformation
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis.DataCorrection\
    import CorrectionByFFT

from FitUtil.WormLikeChain.Python.Code.WLC_Fit import BoundedWlcFit
from FitUtil.FitUtils.Python.FitClasses import GetBoundsDict
from GeneralUtil.python import PlotUtilities as pPlotUtil
from GeneralUtil.python import GenUtilities as pGenUtil
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
from GeneralUtil.python.IgorUtil import SavitskyFilter

def ReadInData(FullName):
    """
    Given a file name, reads in its data as TimeSepForce Objects
    
    Args:
        FullName: full path to the file
    Returns:
        List of TimeSepForce Objects
    """
    mObjs = FEC_Util.ReadInData(FullName)
    return mObjs

def GetCorrectedFECs(DataArray):
    """
    Given a (raw) data, gets the corrected (ie: zero offset, flipped, and 
    interference corrected) data
    
    Args:
         DataArray: output of ReadInData
    Returns:
         List, each element is tuple of <Approach,Retract> which are the
         corresponding parts of the curve
    """
    Corrected = []
    for i,Tmp in enumerate(DataArray):
        CorrectionObj = CorrectionByFFT.CorrectionObject()
        Approach,Retract = CorrectionObj.CorrectApproachAndRetract(Tmp)
        # zero out the forces, pre-process...
        ApproachCorrected,RetractCorrected =\
            FEC_Util.PreProcessApproachAndRetract(Approach,Retract)
        Corrected.append([ApproachCorrected,RetractCorrected])
    return Corrected

def GetWLCFits(CorrectedApproachAndRetracts,NoTriggerDistance,
               MaxContourLength,Ns=30):
    """
    Gets the WLC Fits and associated FirObjects

    Args:
        see GetTransitionForces, except...
        MaxContourLength: maximum contour length to allow
        Ns: number of points on the Lp vs L0 grid (fitting grid)
    Returns:
        List of tuples, each tuple is like <Separation for WLC fit, FitObject>
    """
    # get all the fits to the corrected WLCs
    ToRet = []
    for i,(Approach,Retract) in enumerate(CorrectedApproachAndRetracts):
        # get the zerod and corrected Force extension curve
        WLCFitFEC = FEC_Util.\
            GetRegionForWLCFit(Retract,NoTriggerDistance=NoTriggerDistance)
        # actually fit the WLC
        Bounds = GetBoundsDict(**dict(Lp=[0.3e-9,60e-9],
                                      L0=[0,MaxContourLength],
                                      K0=[1000e-12,1400e-12],
                                      kbT=[0,np.inf]))
        SepNear = WLCFitFEC.Separation
        ForceNear = WLCFitFEC.Force
        Fit = BoundedWlcFit(SepNear,ForceNear,VaryL0=True,VaryLp=True,Ns=Ns,
                            Bounds=Bounds)
        Pred = Fit.Predict(SepNear)
        ToRet.append([SepNear,Fit])
        print("{:d}/{:d}".format(i,len(CorrectedApproachAndRetracts)))
        print(Fit)
    return ToRet

def GetTransitionForces(CorrectedApproachAndRetracts,
                        NoTriggerDistance,
                        ContourLength):
    """
    Gets the transition forces associated with the WLC curves. We assume they
    make it all the way through the first and second WLC portions (ie: 
    two enthaulpic regions

    Args:
        CorrectedApproachAndRetracts: List of tuples, each of which is 
        <Corrected Approach, Corrected Retract> TimeSepForce Objectd
      
        NoTriggerDistance: distance where adhesions may be happening; ignore 
        before
        
        ContourLength: 
    Returns:
        List, each element is (hopefully) the transition region
    """
    ToRet = []
    for i,(Approach,Retract) in enumerate(CorrectedApproachAndRetracts):
        # get the (entire) zerod and corrected retract curve
        NFilterPoints = 10
        Retract =  FEC_Util.\
            GetFECPullingRegionAlreadyFlipped(Retract,
                                              NFilterPoints=NFilterPoints)
        # get the normal points and outliers
        AdhesionArgs = dict(Retract=Retract,
                            NFilterPoints=NFilterPoints,
                            NoTriggerDistance=NoTriggerDistance)
        WlcObj =  FEC_Util.\
                  GetGradientOutliersAndNormalsAfterAdhesion(**AdhesionArgs)
        # get the WLC index object
        Idx = FEC_Util.GetWlcIdxObject(WlcObj,Retract)
        # the transition point is from the end of the first WLC to the start
        # of the secon
        StartIdx = Idx.FirstWLC.end
        EndIdx = Idx.SecondWLC.start
        # get that region
        TransitionRegionSlice =slice(StartIdx,EndIdx,1)
        RetractTx = Idx.TimeSepForceObject
        TransitionForce = RetractTx.Force[TransitionRegionSlice]
        TransitionSeparation = RetractTx.Separation[TransitionRegionSlice]
        ToRet.append(TransitionForce)
    return ToRet

def PlotFits(Corrected,ListOfSepAndFits,TransitionForces,ExpectedContourLength):
    """
    Plots the WLC fits into ./Out/

    Args:
        Corrected: list of tuples, each of which is an approach/corrected
        TimeSepForce Object 
        ListOfSepAndFits: see output of GetWLCFits
        TransitionForces: see output of GetTransitionForces
        ExpectedContourLength: expected contour length, in meters
    """
    # get all of the transition forces
    ExpectedContourLengthNm = ExpectedContourLength * 1e9
    # POST: have everything corrected, fit...
    MaxX = ExpectedContourLengthNm*3
    set_y_lim = lambda :  plt.ylim([-50,170])
    set_x_lim = lambda :  plt.xlim([-10,MaxX])
    LegendOpts = dict(loc='upper right',frameon=True)
    # note: we want to pre-process (convert to sensible units etc) but no
    # need to correct (ie: flip and such)
    FilterSpatialResolution = 0.5e-9
    for i,(ApproachCorrected,RetractCorrected) in enumerate(Corrected):
        SepNear,FitObj = ListOfSepAndFits[i]
        # get the number of filter points needed for whatever spatial resolution
        # we want
        NFilterPoints = 75
        PlotOptions = dict(LegendOpts=LegendOpts,
                           NFilterPoints=NFilterPoints)
        # get the WLC prediction for the region we fit
        SepPred = np.linspace(0,max(SepNear),num=50)
        WLC_Pred = FitObj.Predict(SepPred)
        # convert to plotting units
        ToYUnits = lambda y : y*1e12
        ToXUnits = lambda x: x*1e9
        WLC_Force_pN = ToYUnits(WLC_Pred)
        WLC_Separation_nm = ToXUnits(SepPred)
        # Get the fit parameters
        L0,Lp,_,_ = FitObj.Params()
        fig = pPlotUtil.figure()
        SaveNameIncremental = lambda j : "./Out/FEC{:d}_{:d}.png".format(i,j)
        FEC_Plot.FEC_AlreadySplit(ApproachCorrected,RetractCorrected,
                                  **PlotOptions)
        set_y_lim()
        set_x_lim()
        pPlotUtil.LegendAndSave(fig,SaveNameIncremental(0),**LegendOpts)
        plt.plot(WLC_Separation_nm,
                 WLC_Force_pN,linewidth=1.5,color='g',linestyle='-',
                 label="WLC: L0={:4.1f}nm".format(L0*1e9))
        pPlotUtil.LegendAndSave(fig,SaveNameIncremental(1),**LegendOpts)
        ContourLengthLabel = r"""L$_0$={:4.1f}nm""".\
                             format(ExpectedContourLengthNm)
        plt.axvline(ExpectedContourLengthNm,label=ContourLengthLabel,
                    linewidth=5.0,color='g',linestyle='--')
        plt.axhline(65,label=r'F$_{\rm Overstretch}$=65pN',
                    linewidth=5.0,color='k',linestyle='-')
        plt.axhline(ToYUnits(np.median(TransitionForces[i])),linestyle='--')
        pPlotUtil.legend(**LegendOpts)
        pPlotUtil.savefig(fig,SaveNameIncremental(2),**LegendOpts)
        plt.close(fig)


def ScatterPlot(TransitionForces,ListOfSepAndFits,ExpectedContourLength):
    """
    Makes a scatter plot of the contour length and transition forces

    Args:
        TransitionForces: array, each element the transition region for curve i
        ListOfSepAndFits: array, each element the output of GetWLCFits
    """
    L0Arr = []
    TxArr = []
    for (SepNear,FitObj),TransitionFoces in zip(ListOfSepAndFits,
                                                TransitionForces):
        MedianTx = np.median(TransitionFoces)
        L0,Lp,_,_ = FitObj.Params()
        L0Arr.append(L0)
        TxArr.append(MedianTx)
    # go ahead an throw out ridiculous data from the WLC, where transition
    # normalize the contour length to L0
    L0Arr = np.array(L0Arr)/ExpectedContourLength
    # convert to useful units
    L0Plot = np.array(L0Arr)
    TxPlot =  np.array(TxArr) * 1e12
    fig = pPlotUtil.figure(figsize=(12,12))
    plt.subplot(2,2,1)
    plt.plot(L0Plot,TxPlot,'go',label="Data")
    alpha = 0.3
    ColorForce = 'r'
    ColorLength = 'b'
    plt.axhspan(62,68,color=ColorForce,label=r"$F_{\rm tx}$ $\pm$ 5%",
                alpha=alpha)
    L0BoxMin = 0.9
    L0BoxMax = 1.1
    plt.axvspan(L0BoxMin,L0BoxMax,color=ColorLength,
                label=r"L$_{\rm 0}$ $\pm$ 10%",alpha=alpha)
    fudge = 1.05
    # make the plot boundaries OK
    MaxX = max(L0BoxMax,max(L0Plot))*fudge
    MaxY = max(max(TxPlot),68)*fudge
    plt.xlim([0,MaxX])
    plt.ylim([0,MaxY])
    pPlotUtil.lazyLabel("",r"F$_{\rm overstretch}$ (pN)",
                        "DNA Characterization Histograms ",frameon=True)
    ## now make 1-D histograms of everything
    # subplot of histogram of transition force
    HistOpts = dict(alpha=alpha,linewidth=0)
    plt.subplot(2,2,2)
    TransitionForceBins = np.linspace(0,MaxY)
    plt.hist(TxPlot,bins=TransitionForceBins,orientation="horizontal",
             color=ColorForce,**HistOpts)
    pPlotUtil.lazyLabel("Count","","")
    plt.ylim([0,MaxY])
    plt.subplot(2,2,3)
    ContourBins = np.linspace(0,MaxX)
    plt.hist(L0Plot,bins=ContourBins,color=ColorLength,**HistOpts)
    pPlotUtil.lazyLabel(r"$\frac{L_{\rm WLC}}{L_0}$","Count","")
    plt.xlim([0,MaxX])
    pPlotUtil.savefig(fig,"./Out/ScatterL0vsFTx.png")

def run():
    """
    Runs contour length analysis
    """
    DataDir ="./Data/"
    FullNames = pGenUtil.getAllFiles(DataDir,".pxp")
    DataArray = []
    Force = False
    ForceWLC = False
    ForceTransition = False
    MetersPerBp = 0.338e-9
    # primer locations, plus overhang, plus abasic sites
    Bp = 3520-1607+12+3
    ExpectedContourLength =MetersPerBp*Bp
    NoTriggerDistance = ExpectedContourLength/4
    GridResolution = 200
    # maximum contour length should take into account linkers, tip (~30nm)
    WlcParams = dict(MaxContourLength=ExpectedContourLength*1.1+30e-9,
                     NoTriggerDistance=NoTriggerDistance,
                     Ns=GridResolution)
    # read in all the data
    for i,Name in enumerate(FullNames):
        FileName = Name[len(DataDir):]
        DataArray.extend(pCheckUtil.getCheckpoint("./Cache/{:s}.pkl".\
                                                  format(FileName),
                                                  ReadInData,Force,Name))
    # get all the corrected, force-zeroed objects
    Corrected= pCheckUtil.getCheckpoint("./Cache/Corrected.pkl",
                                        GetCorrectedFECs,
                                        Force,DataArray)
    # Get all the WLC (initial)
    ListOfSepAndFits= pCheckUtil.getCheckpoint("./Cache/WLC.pkl",GetWLCFits,
                                               (ForceWLC),
                                               Corrected,**WlcParams)
    TransitionArgs = dict(NoTriggerDistance=NoTriggerDistance,
                          CorrectedApproachAndRetracts=Corrected,
                          ContourLength=ExpectedContourLength)
    TransitionForces= pCheckUtil.getCheckpoint("./Cache/Transition.pkl",
                                               GetTransitionForces,
                                               ForceTransition,
                                               **TransitionArgs)
    # make a scatter plot of L0 and the overstrectching force. 
    ScatterPlot(TransitionForces,ListOfSepAndFits,
                ExpectedContourLength=ExpectedContourLength)
    # plot the WLC fits
    PlotFits(Corrected,ListOfSepAndFits,TransitionForces,
             ExpectedContourLength=ExpectedContourLength)

if __name__ == "__main__":
    run()
