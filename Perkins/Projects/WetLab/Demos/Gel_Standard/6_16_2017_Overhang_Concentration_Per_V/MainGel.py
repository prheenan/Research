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
from collections import OrderedDict

class DistributionObj:

    def __init__(self,Linear,Circular,Concat):
        self.Linear = Linear
        self.Circular = Circular
        self.Concat = Concat

class OverhangLane:
    """
    Class to keep track of the bands in an overhang lane
    """
    def __init__(self,Linear,Circular,*Concat):
        self.LinearBand=Linear
        self.CircularBand = Circular
        self.Concatemers = Concat
        self.Total = Linear + Circular + sum(Concat)
    def _Norm(self,x):
        return x/self.Total
    @property
    def LinearRelative(self):
        return self._Norm(self.LinearBand)
    @property
    def CircularRelative(self):
        return self._Norm(self.CircularBand)
    @property
    def ConcatemerRelative(self):
        return self._Norm(sum(self.Concatemers))
    def __str__(self):
        return "Lin:{:3.2f},Circ:{:3.2f},Concat:{:3.2f}".\
            format(self.LinearRelative,
                   self.CircularRelative,
                   self.ConcatemerRelative)
    def __repr__(self):
        return str(self)


def GetImageJData(DataDirBase,ext=".xls"):
    """
    Given a base data directory, finds all files with ext in each subdirectory

    Args:
        DataDirBase: base data directory. Each subdirectory has files with 
        extension 'ext'
       
        ext: file extension
    Returns:
        ordered dictionary of <subdir:fullpaths>
    """
    Voltages = OrderedDict()
    for f in sorted(os.listdir(DataDirBase)):
        PossibleSubDir = DataDirBase + f +"/"
        if (os.path.isdir(PossibleSubDir)):
            Files = pGenUtil.getAllFiles(PossibleSubDir,".xls")
            Voltages[f] =Files
    return Voltages

def GetImageJMeasurements(File):
    """
    Returns the in-order values of the intensity column in the ImageJ xls file

    Args:
        File: to read from
    Returns:
        intensity column
    """
    return np.loadtxt(File,skiprows=1,usecols=(1,))

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
            # get this lanes data, reverse so it goes linear,circular,conat
            DataTmp = GetImageJMeasurements(f)[::-1]
            ThisLane = OverhangLane(*DataTmp)
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
                        "DNA Populations as a function of voltage",frameon=True)
    plt.xlim([0,np.max(VoltageFloats)*1.1])
    pPlotUtil.savefig(fig,"BoxPlots.png")
    
if __name__ == "__main__":
    run()
