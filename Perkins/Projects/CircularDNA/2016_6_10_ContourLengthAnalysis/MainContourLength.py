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

from FitUtil.WormLikeChain.Python.Code.WLC_Fit import BoundedWlcFit
from FitUtil.FitUtils.Python.FitClasses import GetBoundsDict
from GeneralUtil.python import PlotUtilities as pPlotUtil
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
import copy
from GeneralUtil.python.IgorUtil import SavitskyFilter


def ReadInData(FullName):
    mObjs = FEC_Util.ReadInData(FullName)
    return mObjs

def run():
    """
    Runs contour length analysis
    """
    OutFile = ""
    Limit = 2
    FullName = "2016_7_5_hot_depositition_" +\
               "water-sealed-better-efficiency-still-non-circular.pxp"
    DataArray = pCheckUtil.getCheckpoint("Tmp.pkl",ReadInData,False,FullName)
    NoTriggerDistance = 100e-9
    for Tmp in DataArray:
        idx = 0
        Corrected,_ = CorrectForcePullByMetaInformation(Tmp)
        Sep = Tmp.Separation
        Tmp = Corrected
        # work with the corrected version
        Approach,Retract = FEC_Util.GetApproachRetract(Tmp)
        EntireRetract = FEC_Util.\
                        GetFECPullingRegion(Retract,MetersAfterTouchoff=None)
        FilterFactor =10
        NFilterPoints = int(np.ceil(EntireRetract.Force.size/FilterFactor))
        FilteredForce = FEC_Util.GetFilteredForce(EntireRetract,NFilterPoints)
        FilteredForceGradient = SavitskyFilter(np.gradient(FilteredForce.Force),
                                               NFilterPoints)
        OnlyPositive = FilteredForceGradient[np.where(FilteredForceGradient>0)]
        q75, q25 = np.percentile(OnlyPositive, [75 ,25])
        iqr = q75-q25
        IsOutlier = lambda x: x > q75 + 1.5 * iqr
        FirstOutlier = np.where(IsOutlier(FilteredForceGradient))[0][0]
        MaxIdx = np.argmax(FilteredForceGradient)
        IdxArr = np.arange(0,FilteredForceGradient.size)
        SeparationRelativeRetract = EntireRetract.Separation
        SeparationRelativeRetract -= SeparationRelativeRetract[0]
        # first worm like chain ends where we past the max no longer an outlier
        Outliers = np.where( ~IsOutlier(FilteredForceGradient) &
                             (IdxArr > FirstOutlier) &
                             (SeparationRelativeRetract > NoTriggerDistance))
        EndOfFirstWLC = Outliers[0][0]
        MetersAfterTouchoff = SeparationRelativeRetract[EndOfFirstWLC]
        NearSurface =  FEC_Util.\
                GetFECPullingRegion(Retract,
                                    MetersAfterTouchoff=MetersAfterTouchoff)
        Bounds = GetBoundsDict(**dict(Lp=[20e-9,60e-9],
                                      L0=[100e-9,700e-9],
                                      K0=[1000e-12,1400e-12],
                                      kbT=[0,np.inf]))
        SepNear = NearSurface.Separation
        ForceNear = NearSurface.Force
        Fit = BoundedWlcFit(SepNear,ForceNear,VaryL0=True,VaryLp=True,Ns=20,
                            Bounds=Bounds)
        Pred = Fit.Predict(SepNear)
        # the fit was to 'NearSurface', which is zeroed. There can be an offset
        # due to hydrodynamic drag on the cantilever (typically <20pN)
        # to find this (in order for the retract-offsetted WLC fit to match),
        # we each
        Appr,Retr = FEC_Util.SplitAndProcess(Tmp)
        # how much of the retract should we use to figure out the zero?
        fraction = 0.05
        N = int(np.ceil(fraction*Retr.Force.size))
        # get the two zeros, and offset the fit by their different (retract
        # should almost certainly be higher)
        ZeroAppr = np.median(Appr.Force[:N])
        ZeroRetr = np.median(Retr.Force[-N:])
        Offset = ZeroRetr - ZeroAppr
        # offset the WLC 
        Pred += Offset
        # plot the data and the prediction
        fig = pPlotUtil.figure()
        FEC_Plot.FEC(Tmp)
        # now plot some meta information. The expected overstretch
        ExpectedOverstretch_pN = 65
        plt.axhline(ExpectedOverstretch_pN,
                    linewidth=3.0,color='k',linestyle="--",label="65pN")
        ToNm = lambda x: x*1e9
        ToPn = lambda x: x*1e12
        # get the contour length (L0) and the persistence length (Lp) in nm
        # and as integers (ie: dont care about 2%
        L0_nm = int(ToNm(Fit.Info.ParamVals.ParamDict["L0"].Value))
        Lp_nm = int(ToNm(Fit.Info.ParamVals.ParamDict["Lp"].Value))
        plt.axvline(ToNm(MetersAfterTouchoff))
        # plot the WLC prediction, label...
        plt.plot(ToNm(SepNear),ToPn(Pred),color='g',linestyle='--',
                 linewidth=5.0,
                 label="WLC (Extensible)\n" +\
                 r"$L_0$={:d}nm, $L_p$={:d}nm".format(L0_nm,Lp_nm))
        pPlotUtil.legend(frameon=True)
        # note: limits are in nm and pN
        MaxY_pN = np.max(ToPn(Pred[np.where(np.isfinite(Pred))]))
        MaxY_pN = max(MaxY_pN,ToPn(np.max(Retr.Force)))
        MinY_pN = -MaxY_pN/5
        plt.ylim([MinY_pN,MaxY_pN])
        plt.xlim([-20,plt.xlim()[-1]])
        Name = "WLC" + Tmp.Meta.Name
        pPlotUtil.savefig(fig,Name + ".png")
        # Read in the pxp (assume each 'name-group' with the same numerical
        # suffix represents a valid wave with a WLC of interest)
    
    

if __name__ == "__main__":
    run()
