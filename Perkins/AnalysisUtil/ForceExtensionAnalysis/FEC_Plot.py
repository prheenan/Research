# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import FEC_Util
import GeneralUtil.python.PlotUtilities as pPlotUtil

import copy


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

def FEC_AlreadySplit(Appr,Retr,
                     XLabel = "Separation (nm)",
                     YLabel = "Force (pN)",
                     ConversionOpts=dict(ConvertX = lambda x: x*1e9,
                                         ConvertY = lambda y: y*1e12),
                     PlotLabelOpts=dict(),
                     PreProcess=False,
                     NFilterPoints=50,
                     LegendOpts=dict(loc='best'),
                     **kwargs):
    """

    Args:
        XLabel: label for x axis
        YLabel: label for y axis
        ConversionOpts: see FEC_Util.SplitAndProcess
        PlotLabelOpts: see arguments after filtering of ApproachRetractCurve
        PreProcess: if true, pre-processes the approach and retract separately
        (ie: to zero and flip the y axis).
        NFilterPoints: see FEC_Util.SplitAndProcess, for Savitsky-golay
        PreProcess: passed to 

    """
    ApprCopy = FEC_Util.UnitConvert(Appr,**ConversionOpts)
    RetrCopy = FEC_Util.UnitConvert(Retr,**ConversionOpts)
    if (PreProcess):
        ApprCopy,RetrCopy = FEC_Util.PreProcessApproachAndRetract(ApprCopy,
                                                                  RetrCopy,
                                                                  **kwargs)
    _ApproachRetractCurve(ApprCopy,RetrCopy,
                          NFilterPoints=NFilterPoints,**PlotLabelOpts)
    pPlotUtil.lazyLabel(XLabel,YLabel,"")
    pPlotUtil.legend(**LegendOpts)
    
def FEC(TimeSepForceObj,
        PreProcessDict=dict(),
        **kwargs):
    """
    Plots a force extension curve. Splits the curve into approach and 
    Retract and pre-processes by default

    Args:
        TimeSepForceObj: 'Raw' TimeSepForce Object
        PreProcessDict: passed directly to FEC_Util.PreProcessFEC
        **kwargs: passed directly to FEC_Plot.FEC_AlreadySplit
    """
    Appr,Retr= FEC_Util.PreProcessFEC(TimeSepForceObj,
                                      NFilterPoints=NFilterPoints,
                                      **PreProcessDict)
    # plot the approach and retract with the appropriate units
    FEC_AlreadySplit(Appr,Retr,**kwargs)
