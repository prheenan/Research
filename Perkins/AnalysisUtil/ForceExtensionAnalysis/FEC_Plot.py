# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import FEC_Util
import GeneralUtil.python.PlotUtilities as pPlotUtil

import copy

def PreProcessFEC(TimeSepForceObject,NFilterPoints=100,
                  ZeroForceFraction=0.2,
                  ZeroSep=True,FlipY=True):
    """
        ZeroForceFraction: if not None, fraction of points near the retract end
        to filter to
        
        ZeroSep: if true, zeros the separation to its minima
        FlipY: if true, multiplies Y (force) by -1 before plotting

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
    return Appr,Retr

def _ApproachRetractCurve(Appr,Retr,NFilterPoints=100,
                          ApproachLabel="Approach",
                          RetractLabel="Retract"):
    """
    Most of the brains for the approach/retract curve. does *not* show anything

    Args:
        TimeSepForceObject: what we are plotting
        NFilterPoints: how many points to filter down
        ApproachLabel: label to put on the approach
        RetractLabel: label to put on the retract
    """
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
    """
    Converts the 'X' and 'Y' using the specified units and properties 
    of the object passed in 

    Args:
        TimeSepForceObj : see ApproachRetractCurve
        ConvertX: method to convert the X values into whatever units we want
        ConvertY: metohod to convery the Y values into whatever units we want
        GetX: gets the x values (assumed separation for plotting XXX TODO)
        GetY: gets the y values (assumed force for plotting XXX TODO)
    Returns: 
        deep *copy* of original object in the specified units
    """
    ObjCopy = copy.deepcopy(TimeSepForceObj)
    ObjCopy.Force = ConvertY(GetY(ObjCopy))
    ObjCopy.Separation = ConvertX(GetX(ObjCopy))
    return ObjCopy
    
def FEC(TimeSepForceObj,
        XLabel="Separation (nm)",
        YLabel="Force (pN)",
        ConversionOpts=dict(),
        PlotLabelOpts=dict(),
        NFilterPoints=50,
        **kwargs):
    
    # convert the x and y to sensible units
    ObjCopy = UnitConvert(TimeSepForceObj,**ConversionOpts)
    # pre-process (to, for example, flip the axes and zero everything out
    Appr,Retr = PreProcessFEC(ObjCopy,NFilterPoints=NFilterPoints,**kwargs)
    # actually plot, with the filtered versions
    _ApproachRetractCurve(Appr,Retr,NFilterPoints=NFilterPoints,**PlotLabelOpts)
    pPlotUtil.lazyLabel(XLabel,YLabel,"")
    pPlotUtil.legend()
