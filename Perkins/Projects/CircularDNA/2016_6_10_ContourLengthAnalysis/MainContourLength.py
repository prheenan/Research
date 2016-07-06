# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot

from FitUtil.WormLikeChain.Python.Code.WLC_Fit import BoundedWlcFit
from FitUtil.FitUtils.Python.FitClasses import GetBoundsDict
from GeneralUtil.python import PlotUtilities as pPlotUtil


def run():
    """
    Runs contour length analysis
    """
    OutFile = ""
    Limit = 2
    "ContourLengthExperiment.pxp"
    FullName = "ContourLengthExperiment.pxp"
    mObjs = FEC_Util.ReadInData(FullName)
    # get where the surface of this object is
    Tmp = mObjs[0]
    Approach,Retract = FEC_Util.GetApproachRetract(Tmp)
    NearSurface =  FEC_Util.GetFECPullingRegion(Retract,
                                                MetersAfterTouchoff=640e-9)
    Bounds = GetBoundsDict(**dict(Lp=[35e-9,60e-9],
                                  L0=[600e-9,700e-9],
                                  K0=[1000e-12,1400e-12],
                                  kbT=[0,np.inf]))
    SepNear = NearSurface.Separation
    ForceNear = NearSurface.Force
    Fit = BoundedWlcFit(SepNear,ForceNear,VaryL0=True,VaryLp=True,Ns=20,
                        Bounds=Bounds)
    Pred = Fit.Predict(SepNear)
    # offset the WLC, in pN
    Pred += 10e-12
    # plot the data and the prediction
    fig = plt.figure()
    FEC_Plot.FEC(Tmp,NFilterPoints=40)
    plt.axhline(65,linewidth=3.0,color='k',linestyle="--",label="65pN")
    ToNm = lambda x: x*1e9
    ToPn = lambda x: x*1e12
    plt.plot(ToNm(SepNear),ToPn(Pred),color='g',linestyle='--',linewidth=5.0,
             label="Extensible WLC")
    pPlotUtil.legend(frameon=True)
    # note: limits are in nm and pN
    plt.ylim([-30,150])
    plt.xlim([-10,plt.xlim()[-1]])
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
