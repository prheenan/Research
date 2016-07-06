# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import FEC_Util
import GeneralUtil.python.PlotUtilities as pPlotUtil

import copy

def _ApproachRetractCurve(TimeSepForceObject,NFilterPoints=100,
                          ZeroForceFraction=0.2,
                          ZeroSep=True,FlipY=True,
                          ApproachLabel="Approach",
                          RetractLabel="Retract"):
    """
    Most of the brains for the approach/retract curve. does *not* show anything

    Args:
        TimeSepForceObject: what we are plotting
        NFilterPoints: how many points to filter down
        ZeroForceFraction: if not None, fraction of points near the retract end
        to filter to
        
        ZeroSep: if true, zeros the separation to its minima
        FlipY: if true, multiplies Y (force) by -1 before plotting
        ApproachLabel: label to put on the approach
        RetractLabel: label to put on the retract
    """
    Appr,Retr = FEC_Util.GetApproachRetract(TimeSepForceObject)
    if (ZeroForceFraction is not None):
        # then we need to offset the force
        # XXX assume offset is the same for both
        _,ZeroForceRetr = FEC_Util.\
                          GetSurfaceIndexAndForce(Retr,
                                                  Fraction=ZeroForceFraction,
                                                  FilterPoints=NFilterPoints,
                                                  ZeroAtStart=False)
        _,ZeroForceAppr = FEC_Util.\
                          GetSurfaceIndexAndForce(Appr,
                                                  Fraction=ZeroForceFraction,
                                                  FilterPoints=NFilterPoints,
                                                  ZeroAtStart=True)
        Appr.Force += ZeroForceAppr
        # Do the same for retract
        Retr.Force += ZeroForceAppr
    if (ZeroSep):
        MinSep = np.min(TimeSepForceObject.Separation)
        Appr.Separation -= MinSep
        Retr.Separation -= MinSep
    if (FlipY):
        Appr.Force *= -1
        Retr.Force *= -1
    # plot the separation and force, with their filtered counterparts
    ApprFiltered = FEC_Util.GetFilteredForce(Appr,NFilterPoints)
    RetrFiltered = FEC_Util.GetFilteredForce(Retr,NFilterPoints)
    plt.plot(Appr.Separation,Appr.Force,color='r',alpha=0.3)
    plt.plot(Appr.Separation,ApprFiltered.Force,color='r',label=ApproachLabel)
    plt.plot(Retr.Separation,Retr.Force,color='b',alpha=0.3)
    plt.plot(Retr.Separation,RetrFiltered.Force,color='b',label=RetractLabel)

def UnitConvert(TimeSepForceObj,
                ConvertX=lambda x : x*1e9,
                ConvertY=lambda y : y*1e12,
                GetX = lambda x : x.Separation,
                GetY = lambda x : x.Force):
    ObjCopy = copy.deepcopy(TimeSepForceObj)
    ObjCopy.Force = ConvertY(GetY(ObjCopy))
    ObjCopy.Separation = ConvertX(GetX(ObjCopy))
    return ObjCopy

    
def FEC(TimeSepForceObj,
        XLabel="Separation (nm)",
        YLabel="Force (pN)",
        ConversionOpts=dict(),
        **kwargs):
    # convert the x and y to sensible units
    ObjCopy = UnitConvert(TimeSepForceObj,**ConversionOpts)
    _ApproachRetractCurve(ObjCopy,**kwargs)
    pPlotUtil.lazyLabel(XLabel,YLabel,"")
    pPlotUtil.legend()
