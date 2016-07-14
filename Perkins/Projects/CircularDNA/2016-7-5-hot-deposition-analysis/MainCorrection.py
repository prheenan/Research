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
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
from GeneralUtil.python.IgorUtil import SavitskyFilter


def ReadInData(FullName):
    mObjs = FEC_Util.ReadInData(FullName)
    return mObjs

def GetCorrectedFECs(DataArray):
    Corrected = []
    for i,Tmp in enumerate(DataArray):
        CorrectionObj = CorrectionByFFT.CorrectionObject()
        ApproachCorrected,RetractCorrected = \
                CorrectionObj.CorrectApproachAndRetract(Tmp)
        Corrected.append([ApproachCorrected,RetractCorrected])
    return Corrected

def GetWLCFits(CorrectedApproachAndRetracts):
    # get all the fits to the corrected WLCs
    ToRet = []
    for Approach,Retract in CorrectedApproachAndRetracts:
        _,RetractProcessed = FEC_Util.PreProcessApproachAndRetract(Approach,
                                                                   Retract)
        # get the zerod and corrected Force extension curve
        WLCFitFEC = FEC_Util.GetRegionForWLCFit(RetractProcessed)
        # actually fit the WLC
        Bounds = GetBoundsDict(**dict(Lp=[30e-9,55e-9],
                                      L0=[150e-9,700e-9],
                                      K0=[1000e-12,1400e-12],
                                      kbT=[0,np.inf]))
        SepNear = WLCFitFEC.Separation
        ForceNear = WLCFitFEC.Force
        Fit = BoundedWlcFit(SepNear,ForceNear,VaryL0=True,VaryLp=True,Ns=40,
                            Bounds=Bounds)
        Pred = Fit.Predict(SepNear)
        ToRet.append([SepNear,Fit])
    return ToRet

def run():
    """
    Runs contour length analysis
    """
    OutFile = ""
    Limit = 2
    FullNames = ["2016_7_10_1ng_ul_50C_4hour_depo_circ_dna_Strept_tip_I.pxp"]
    DataArray = []
    # read in all the data
    for i,Name in enumerate(FullNames):
        DataArray.extend(pCheckUtil.getCheckpoint("Tmp{:d}.pkl".format(i),
                                                  ReadInData,False,Name))
    # get all the corrected objects
    Corrected= pCheckUtil.getCheckpoint("Corrected.pkl",GetCorrectedFECs,
                                        False,DataArray)
    # Get all the WLCs
    ListOfSepAndFits= pCheckUtil.getCheckpoint("WLC.pkl",GetWLCFits,
                                               False,Corrected)

    # POST: have everything corrected, fit...
    set_y_lim = lambda :  plt.ylim([-50,200])
    set_x_lim = lambda :  plt.xlim([-10,1300])
    LegendOpts = dict(loc='upper left',frameon=True)
    # note: we want to pre-process (convert to sensible units etc) but no
    # need to correct (ie: flip and such)
    PlotOptions = dict(FlipY=True,
                       PreProcess=True,
                       LegendOpts=LegendOpts)
    for i,(ApproachCorrected,RetractCorrected) in enumerate(Corrected):
        SepNear,FitObj = ListOfSepAndFits[i]
        # get the WLC prediction
        WLC_Pred = FitObj.Predict(SepNear)
        # convert to plotting units
        WLC_Force_pN = WLC_Pred  * 1e12
        WLC_Separation_nm = SepNear * 1e9
        fig = pPlotUtil.figure()
        FEC_Plot.FEC_AlreadySplit(ApproachCorrected,RetractCorrected,
                                  **PlotOptions)
        plt.plot(WLC_Separation_nm,
                 WLC_Force_pN,linewidth=3,color='g',linestyle='--',
                 label="WLC PRediction")
        plt.axhline(65,label="65pN",linewidth=5.0,color='g',linestyle='--')
        set_y_lim()
        set_x_lim()
        pPlotUtil.legend(**LegendOpts)
        pPlotUtil.savefig(fig,"./tmp{:d}.png".format(i))

if __name__ == "__main__":
    run()
