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

def FEC(TimeSepForceObj,
        XLabel="Separation (nm)",
        YLabel="Force (pN)",
        ConversionOpts=dict(ConvertX = lambda x: x*1e9,
                            ConvertY = lambda y: y*1e12),
        PlotLabelOpts=dict(),
        NFilterPoints=50,
        **kwargs):
    """
    Plots a force extension curve. Very flexible

    Args:
        TimeSepForceObj: 'Raw' TimeSepForce Object
        XLabel: label for x axis
        YLabel: label for y axis
        ConversionOpts: see FEC_Util.SplitAndProcess
        PlotLabelOpts: see arguments after filtering of ApproachRetractCurve
        NFilterPoints: see FEC_Util.SplitAndProcess, for Savitsky-golay
        **kwargs: passed directly to FEC_Util.SplitAndProcess
    """
    Appr,Retr= FEC_Util.SplitAndProcess(TimeSepForceObj,
                                        ConversionOpts=ConversionOpts,
                                        NFilterPoints=NFilterPoints,**kwargs)
    # actually plot, with the filtered versions
    _ApproachRetractCurve(Appr,Retr,NFilterPoints=NFilterPoints,**PlotLabelOpts)
    pPlotUtil.lazyLabel(XLabel,YLabel,"")
    pPlotUtil.legend()
