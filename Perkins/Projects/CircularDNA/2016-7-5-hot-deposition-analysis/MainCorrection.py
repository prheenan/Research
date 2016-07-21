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
    mObjs = FEC_Util.ReadInData(FullName)
    return mObjs

def GetCorrectedFECs(DataArray):
    Corrected = []
    for i,Tmp in enumerate(DataArray):
        CorrectionObj = CorrectionByFFT.CorrectionObject()
        Approach,Retract = CorrectionObj.CorrectApproachAndRetract(Tmp)
        # zero out the forces, pre-process...
        ApproachCorrected,RetractCorrected =\
            FEC_Util.PreProcessApproachAndRetract(Approach,Retract)
        Corrected.append([ApproachCorrected,RetractCorrected])
    return Corrected

def GetWLCFits(CorrectedApproachAndRetracts):
    # get all the fits to the corrected WLCs
    ToRet = []
    for i,(Approach,Retract) in enumerate(CorrectedApproachAndRetracts):
        # get the zerod and corrected Force extension curve
        WLCFitFEC = FEC_Util.GetRegionForWLCFit(Retract)
        # actually fit the WLC
        Bounds = GetBoundsDict(**dict(Lp=[0.3e-9,80e-9],
                                      L0=[120e-9,700e-9],
                                      K0=[1000e-12,1400e-12],
                                      kbT=[0,np.inf]))
        SepNear = WLCFitFEC.Separation
        ForceNear = WLCFitFEC.Force
        Fit = BoundedWlcFit(SepNear,ForceNear,VaryL0=True,VaryLp=True,Ns=60,
                            Bounds=Bounds)
        Pred = Fit.Predict(SepNear)
        ToRet.append([SepNear,Fit])
        print("{:d}/{:d}".format(i,len(CorrectedApproachAndRetracts)))
        print(Fit)
    return ToRet

def GetTransitionForces(CorrectedApproachAndRetracts):
    ToRet = []
    for i,(Approach,Retract) in enumerate(CorrectedApproachAndRetracts):
        # get the (entire) zerod and corrected retract curve
        RetractPull =  FEC_Util.GetFECPullingRegionAlreadyFlipped(Retract)
        # get the normal points and outliers
        NoAdhesionMask,Outliers,Normal =  FEC_Util.\
                GetGradientOutliersAndNormalsAfterAdhesion(RetractPull,
                                                           NFilterFraction=0.05)
        # get the WLC index object
        Idx = FEC_Util.GetWlcIdxObject(NoAdhesionMask,Outliers,Normal,
                                       RetractPull)
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

def PlotFits(Corrected,ListOfSepAndFits,TransitionForces):
    # get all of the transition forces 
    # POST: have everything corrected, fit...
    set_y_lim = lambda :  plt.ylim([-50,150])
    set_x_lim = lambda :  plt.xlim([-10,1500])
    LegendOpts = dict(loc='upper left',frameon=True)
    # note: we want to pre-process (convert to sensible units etc) but no
    # need to correct (ie: flip and such)
    PlotOptions = dict(LegendOpts=LegendOpts)
    for i,(ApproachCorrected,RetractCorrected) in enumerate(Corrected):
        SepNear,FitObj = ListOfSepAndFits[i]
        # get the WLC prediction
        WLC_Pred = FitObj.Predict(SepNear)
        # convert to plotting units
        ToYUnits = lambda y : y*1e12
        ToXUnits = lambda x: x*1e9
        WLC_Force_pN = ToYUnits(WLC_Pred)
        WLC_Separation_nm = ToXUnits(SepNear)
        # Get the fit parameters
        L0,Lp,_,_ = FitObj.Params()
        fig = pPlotUtil.figure()
        FEC_Plot.FEC_AlreadySplit(ApproachCorrected,RetractCorrected,
                                  **PlotOptions)
        plt.plot(WLC_Separation_nm,
                 WLC_Force_pN,linewidth=3,color='g',linestyle='--',
                 label="WLC: L0={:.1f}nm".format(L0*1e9))
        plt.axvline(650,label=r'L$_{\rm Contour}$=650nm',
                    linewidth=5.0,color='g',linestyle='--')
        
        plt.axhline(65,label=r'F$_{\rm Overstretch}$=65pN',
                    linewidth=5.0,color='k',linestyle='-')
        plt.axhline(ToYUnits(np.median(TransitionForces[i])),linestyle='--')
        set_y_lim()
        set_x_lim()
        pPlotUtil.legend(**LegendOpts)
        pPlotUtil.savefig(fig,"./Out/FEC{:d}.png".format(i))

def run():
    """
    Runs contour length analysis
    """
    DataDir ="./Data/"
    FullNames = pGenUtil.getAllFiles(DataDir,".pxp")
    DataArray = []
    Force = False
    ForceWLC = False
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
                                               ForceWLC,Corrected)
    TransitionForces= pCheckUtil.getCheckpoint("./Cache/Transition.pkl",
                                               GetTransitionForces,
                                               ForceWLC,Corrected)
    PlotFits(Corrected,ListOfSepAndFits,TransitionForces)

if __name__ == "__main__":
    run()
