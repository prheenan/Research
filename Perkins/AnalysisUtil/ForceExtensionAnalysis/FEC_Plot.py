# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import FEC_Util
import GeneralUtil.python.PlotUtilities as pPlotUtil

import copy

def _ApproachRetractCurve(TimeSepForceObject,FilterN=100,ZeroForceFraction=0.1,
                          ZeroSep=True,FlipY=True,
                          ApproachLabel="Approach",
                          RetractLabel="Retract"):
    Appr,Retr = FEC_Util.GetApproachRetract(TimeSepForceObject)
    if (ZeroForceFraction is not None):
        # then we need to offset the force
        # XXX assume offset is the same for both
        _,ZeroForce = FEC_Util.GetSurfaceIndex(TimeSepForceObject,
                                               Fraction=ZeroForceFraction,
                                               FilterPoints=FilterN)
        Appr.Force += ZeroForce
        Retr.Force += ZeroForce
    if (ZeroSep):
        MinSep = np.min(TimeSepForceObject.Separation)
        Appr.Separation -= MinSep
        Retr.Separation -= MinSep
    if (FlipY):
        Appr.Force *= -1
        Retr.Force *= -1
    # plot the separation and force, with their filtered counterparts
    ApprFiltered = FEC_Util.GetFilteredForce(Appr,FilterN)
    RetrFiltered = FEC_Util.GetFilteredForce(Retr,FilterN)
    plt.plot(Appr.Separation,Appr.Force,color='r',alpha=0.3)
    plt.plot(Appr.Separation,ApprFiltered.Force,color='r',label=ApproachLabel)
    plt.plot(Retr.Separation,Retr.Force,color='b',alpha=0.3)
    plt.plot(Retr.Separation,RetrFiltered.Force,color='b',label=RetractLabel)

def FEC(TimeSepForceObj,
        ConvertX=lambda x : x*1e9,
        ConvertY=lambda y : y*1e12,
        XLabel="Separation (nm)",
        YLabel="Force (pN)",
        **kwargs):
    # convert the x and y
    ObjCopy = copy.deepcopy(TimeSepForceObj)
    ObjCopy.Force = ConvertY(ObjCopy.Force)
    ObjCopy.Separation = ConvertX(ObjCopy.Separation)
    _ApproachRetractCurve(ObjCopy,**kwargs)
    pPlotUtil.lazyLabel(XLabel,YLabel,"")
    pPlotUtil.legend()
