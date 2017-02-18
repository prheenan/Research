# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import FEC_Util
import GeneralUtil.python.PlotUtilities as PlotUtilities

import copy

def_conversion_opts =dict(ConvertX = lambda x: x*1e9,
                          ConvertY = lambda y: y*1e12)

def _ApproachRetractCurve(Appr,Retr,NFilterPoints=100,
                          x_func = lambda x: x.Separation,
                          y_func = lambda y: y.Force, 
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
    plt.plot(x_func(Appr),y_func(Appr),color='r',alpha=0.3)
    plt.plot(x_func(ApprFiltered),y_func(ApprFiltered),color='r',
             label=ApproachLabel)
    plt.plot(x_func(Retr),y_func(Retr),color='b',alpha=0.3)
    plt.plot(x_func(RetrFiltered),y_func(RetrFiltered),color='b',
             label=RetractLabel)

def FEC_AlreadySplit(Appr,Retr,
                     XLabel = "Separation (nm)",
                     YLabel = "Force (pN)",
                     ConversionOpts=def_conversion_opts,
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
    PlotUtilities.lazyLabel(XLabel,YLabel,"")
    PlotUtilities.legend(**LegendOpts)
    
def FEC(TimeSepForceObj,NFilterPoints=50,
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
    
def heat_map_fec(time_sep_force_objects,num_bins=(100,100),
                 ConversionOpts=def_conversion_opts):
    """
    Plots a force extension curve. Splits the curve into approach and 
    Retract and pre-processes by default

    Args:
        TimeSepForceObj: 'Raw' TimeSepForce Object
        PreProcessDict: passed directly to FEC_Util.PreProcessFEC
        **kwargs: passed directly to FEC_Plot.FEC_AlreadySplit
    """                 
    # convert everything...
    objs = [FEC_Util.UnitConvert(r,**ConversionOpts) 
            for r in time_sep_force_objects]
    filtered_data = [(retr.Separation,retr.Force) for retr in objs]
    separations = np.concatenate([r[0] for r in filtered_data])
    forces = np.concatenate([r[1] for r in filtered_data])
    # make a heat map, essentially
    counts, xedges, yedges, Image = plt.hist2d(separations, forces,
                                               bins=num_bins,cmap='afmhot')
    PlotUtilities.lazyLabel("Separation [nm]",
                            "Force [pN]",
                            "Two-Dimensional Force-Separation Histogram")
    cbar = plt.colorbar()
    label = '# of points in (Force,Separation) Bin'
    cbar.set_label(label,labelpad=10,rotation=270) 
