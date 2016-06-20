# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../..")
import GeneralUtil.python.PlotUtilities as pPlotUtil
import GeneralUtil.python.GenUtilities as pGenUtil
import os
from Research.Perkins.AnalysisUtil.Gels.ImageJUtil import \
    GetImageJData,GetImageJMeasurements,ReadFileToOverhangObj
from collections import OrderedDict

class DistributionObj:
    def __init__(self,Linear,Circular,Concat):
        self.Linear = Linear
        self.Circular = Circular
        self.Concat = Concat

    
def ConvertToOverhangObjects(Voltages):
    """
    Given a dictionary of voltage:Filenames pairs (filenames are the lanes),
    converts each to overhang ovjects, returning the corresponding dictionary
    
    Args:
       Voltages: see return of GetImageJData
    Returns:
       Dictionary like <voltages:OverhangLane Objects>
    """
    ObjectsByVoltage = OrderedDict()
    for Volt,Files in Voltages.items():
        TmpObjs = []
        # go through each file (lane)
        for f in Files:
            ThisLane = ReadFileToOverhangObj(f)
            TmpObjs.append(ThisLane)
        # record all the lanes for this object
        ObjectsByVoltage[Volt] = TmpObjs
    return ObjectsByVoltage

def run():
    """
    Shows how the concatemeters etc change with voltage
    """
    DataDirBase = "./DataByVoltage/"
    # get all the files to use
    Voltages = GetImageJData(DataDirBase,ext=".xls")
    # convert the files to data
    OverhangObjsByVoltage = ConvertToOverhangObjects(Voltages)
    # 'flatten' so we have a distribution of linear, circular, and concatemeter
    # intensities at each voltage
    Distributions = OrderedDict()
    for Volt,AllLanes in OverhangObjsByVoltage.items():
        # get every linear, circular, and concatemer band
        AllLin = [l.LinearRelative for l in AllLanes]
        AllCirc = [l.CircularRelative for l in AllLanes]
        AllConcat = [l.ConcatemerRelative for l in AllLanes]
        Distributions[Volt] = DistributionObj(AllLin,AllCirc,AllConcat)
    # plot the distributions at each voltage
    Strings = Distributions.keys()
    VoltageFloats = [float(v) for v in Strings]
    SortIdx = np.argsort(VoltageFloats)
    VoltageFloats = [VoltageFloats[i] for i in SortIdx]
    # get the voltage (X) values)
    Voltages = [Strings[i] for i in SortIdx]
    # get the flattened lists we want
    GetFlat = lambda f: [f(Distributions[k]) for k in Voltages]
    FlatLin = GetFlat(lambda x: x.Linear)
    FlatCirc = GetFlat(lambda x: x.Circular)
    FlatConcat = GetFlat(lambda x: x.Concat)
    CommonStyle = dict(markersize=10)
    LinProps = dict(marker='x',
                    color='g',linewidth=2,label="Linear",**CommonStyle)
    CirProps = dict(marker='o',
                    color='r',linestyle='--',label="Circular",**CommonStyle)
    ConcatProps = dict(color='k',marker='*',
                       linestyle='--',label="Dimers+",**CommonStyle)
    fig = pPlotUtil.figure()
    Mean = lambda dist: [np.mean(v) for v in dist]
    Stdev = lambda dist: [np.std(v) for v in dist]
    PlotDist = lambda dist,props: plt.errorbar(x=VoltageFloats,
                                               y=Mean(dist),
                                               yerr=Stdev(dist),
                                               **props)
    PlotDist(FlatLin,LinProps)
    PlotDist(FlatCirc,CirProps)
    PlotDist(FlatConcat,ConcatProps)
    pPlotUtil.lazyLabel("Voltage (V)","Intensity Fraction",
                        "Circular DNA fraction saturates at low voltage",
                        frameon=True)
    plt.xlim([0,np.max(VoltageFloats)*1.1])
    pPlotUtil.savefig(fig,"BoxPlots.png")
    
if __name__ == "__main__":
    run()
