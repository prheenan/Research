# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


sys.path.append("../../../../..")
# XXX move to utility?

from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
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
                               Extension=o.Separation,
                               Force=o.Force,
                               SpringConstant=o.SpringConstant,
                               Velocity=o.Velocity)
            for o in TimeSepForceObjects]
    InverseWeierstrass.SetAllWorkOfObjects(Objs)
    return Objs


def run():
    Base = "./Data/"
    Limit=50
    FullName = Base + "2016-6-5-micah-1ppm-biolever-long-strept-saved-data.pxp"
    mObjs = pCheckUtil.getCheckpoint(Base + "cache.pkl",FEC_Util.ReadInData,
                                     False,FullName,Limit)
    ApproachList,RetractList = FEC_Util.BreakUpIntoApproachAndRetract(mObjs)
    PastZeroExt = 5e-9
    FilterPoints = 30
    # get just after the touchoff
    Touchoff = FEC_Util.GetAroundTouchoff(RetractList,FilterPoints=30,
                                          MetersAfterTouchoff=PastZeroExt)
    # get the IWT transform objects
    IwtObjects = ToIWTObjects(Touchoff)
    # get the IWT
    LandscapeObj = InverseWeierstrass.FreeEnergyAtZeroForce(IwtObjects,
                                                            NumBins=25)
    fig = plt.figure()
    for Retract,Touch in zip(RetractList,Touchoff):
        ForZ = Retract
        toNano = lambda x : x * 1e9
        toPn = lambda x: x * 1e12
        RetractZ,_ = Retract.ZsnsrAndDeflV
        RetractZ -= np.min(RetractZ)
        plt.subplot(4,1,1)
        plt.plot(toNano(RetractZ),toPn(Retract.Force),alpha=0.3)
        plt.subplot(4,1,2)
        TouchSep = Touch.Separation
        plt.plot(toNano(TouchSep),toPn(Touch.Force),alpha=0.3)
        plt.subplot(4,1,3)
        NanoExt =toNano(LandscapeObj.Extensions)
        FreeEnergyEq = LandscapeObj.EnergyLandscape
        plt.plot(NanoExt,FreeEnergyEq * LandscapeObj.Beta)
        plt.subplot(4,1,4)
        FOneHalf = 10e-12
        TiltedEnergy = (FreeEnergyEq - LandscapeObj.Extensions * FOneHalf)
        TiltedEnergy -= min(TiltedEnergy[np.isfinite(TiltedEnergy)])
        plt.plot(NanoExt,TiltedEnergy * LandscapeObj.Beta)
        plt.ylim([0,20])
    fig.savefig("out.png")

if __name__ == "__main__":
    run()
