# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../")
import GeneralUtil.python.PlotUtilities as pPlotUtil
import GeneralUtil.python.GenUtilities as pGenUtil
import matplotlib.gridspec as gridspec
import copy

class Trace:
    def __init__(self,Time,Donor,Acceptor,FRET,State,FileName):
        self.Time = Time
        self.Donor = Donor
        self.Acceptor = Acceptor
        self.FRET = FRET
        self.HMM = State
        self.FileName = FileName

def ReadTraces(BaseFolder):
    """
    See: ReadKCl, except with NaCl data instead
    """
    files = sorted(pGenUtil.getAllFiles(BaseFolder,ext=".dat"))
    data = []
    for f in files:
        TmpDat = np.loadtxt(f)
        Columns = TmpDat.shape[1]
        Args = [TmpDat[:,i] for i in range(Columns)]
        TmpColumn = Trace(*Args,FileName=(BaseFolder+f))
        data.append(TmpColumn)
    return data

def PlotSingleTrace(tmp,OutputName):
    fig = pPlotUtil.figure(figsize=(8,12))
    CommonStyle = dict(linewidth=3.0,alpha=0.7)
    # style for lazyLabel
    LazyStyle = dict(loc="upper left",frameon=True)
    ZeroedTime = tmp.Time - tmp.Time[0]
    # function to fudge the y limits to fit the labels
    fudgeFunc = lambda : plt.ylim([-100,1200])
    plt.subplot(2,1,1)
    plt.plot(ZeroedTime,tmp.Donor,"g",label="Donor",**CommonStyle)
    plt.plot(ZeroedTime,tmp.Acceptor,"r",label="Acceptor",**CommonStyle)
    pPlotUtil.lazyLabel("","Channel Intensity","",**LazyStyle)
    fudgeFunc()
    plt.subplot(2,1,2)
    plt.plot(ZeroedTime,tmp.HMM,'b',label="HMM Prediction",alpha=1,
             linewidth=1)
    plt.plot(ZeroedTime,tmp.FRET,'k',label="Raw FRET",alpha=0.3,linewidth=2)
    plt.ylim([-0.1,1.4])
    pPlotUtil.lazyLabel("Time (s)","FRET","",**LazyStyle)
    pPlotUtil.savefig(fig,OutputName)


def run():
    """
    Generates FRET histograms from NaCl and KCl data
    """
    data = ReadTraces("./Data/")
    for i,tmp in enumerate(data):
        FileName = pGenUtil.getFileFromPath(tmp.FileName)
        PlotSingleTrace(tmp,"{:d}Out_".format(i) + FileName +  ".png")

if __name__ == "__main__":
    run()
