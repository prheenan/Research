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
import copy


def run():
    """
    Runs contour length analysis
    """
    OutFile = ""
    Limit = 2
    FullName = "ContourLengthExperiment.pxp"
    mObjs = FEC_Util.ReadInData(FullName)
    # get where the surface of this object is
    Tmp = mObjs[0]
    Corrected, CorrectionInfo = CorrectForcePullByMetaInformation(Tmp)
    # work with the corrected version
    Tmp = Corrected
    Approach,Retract = FEC_Util.GetApproachRetract(Tmp)
    NearSurface =  FEC_Util.GetFECPullingRegion(Retract,
                                                MetersAfterTouchoff=640e-9)
    Bounds = GetBoundsDict(**dict(Lp=[35e-9,60e-9],
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
    # to find this, we each
    Appr,Retr = FEC_Util.SplitAndProcess(Tmp)
    # how much of the retract should we use to figure out the zero?
    fraction = 0.05
    N = int(np.ceil(fraction*Retr.Force.size))
    # get the two zeros, and offset the fit by their different (retract
    # should almost certainly be higher)
    ZeroAppr = np.median(Appr.Force[:N])
    ZeroRetr = np.median(Retr.Force[-N:])
    Offset = ZeroRetr - ZeroAppr
    # offset the WLC, in pN
    Pred += Offset
    # plot the data and the prediction
    fig = plt.figure()
    FEC_Plot.FEC(Tmp,NFilterPoints=40)
    ExpectedOverstretch_pN = 65
    plt.axhline(ExpectedOverstretch_pN,
                linewidth=3.0,color='k',linestyle="--",label="65pN")
    ToNm = lambda x: x*1e9
    ToPn = lambda x: x*1e12
    L0 = int(ToNm(Fit.Info.ParamVals.ParamDict["L0"].Value))
    Lp = int(ToNm(Fit.Info.ParamVals.ParamDict["Lp"].Value))
    plt.plot(ToNm(SepNear),ToPn(Pred),color='g',linestyle='--',linewidth=5.0,
             label="WLC (Extensible)\n" +\
             r"$L_0$={:d}nm, $L_p$={:d}nm".format(L0,Lp))
    pPlotUtil.legend(frameon=True,loc='upper left')
    # note: limits are in nm and pN
    plt.ylim([-30,175])
    plt.xlim([-20,plt.xlim()[-1]])
    pPlotUtil.savefig(fig,"WLC_Tmp.png")
    # Read in the pxp (assume each 'name-group' with the same numerical
    # suffix represents a valid wave with a WLC of interest)

    # get just 'Limit' of the waves

    # get just the retract, normalize the force and separation to zeros.

    # fit an extensible worm-like chain to each, using grid search at a
    # 0.1nm scale for contour length and 0.001nm scale for persistence length
    # within 10% of expected (600-700nm,40-50 for persistence)
    # save this course-grained information

    # plot a histogram of the results

    
    
    

if __name__ == "__main__":
    run()
