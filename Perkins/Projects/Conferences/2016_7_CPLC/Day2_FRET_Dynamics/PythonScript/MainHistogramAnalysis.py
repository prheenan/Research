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

def ReadNaCl(BaseFolder):
    files = sorted(pGenUtil.getAllFiles(BaseFolder,ext=".dat"))
    data = []
    for f in files:
        # only care about third column ('C')
        TmpColumn = np.loadtxt(f,skiprows=52)[:,2]
        data.append(TmpColumn)
    return data

def ReadKCl(BaseFolder):
    files = sorted(pGenUtil.getAllFiles(BaseFolder,ext=".csv"))
    data = []
    for f in files:
        data.append(np.loadtxt(f))
    return data

def AxesDefault(ylim,title,ylabel,UseYTicks):
    ax = plt.gca()
    plt.setp(ax.get_xticklabels(),visible=False)
    if (not UseYTicks):
        plt.setp(ax.get_yticklabels(),visible=False)
    plt.ylim(ylim)
    pPlotUtil.lazyLabel("",ylabel,title,frameon=True)
   
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    data = ReadNaCl("./Data/NaCl")
    dataKcl = ReadKCl("./Data/KCl")
    # assume number of histograms is the same
    NumHists = len(data)
    fig = pPlotUtil.figure(figsize=(12,16))
    CommonStyle = dict(alpha=0.3,linewidth=0)
    StyleDicts = [dict(color='m',label="0mM NaCl",**CommonStyle),
                  dict(color='r',label="1mM NaCl",**CommonStyle),
                  dict(color='k',label="10mM NaCl",**CommonStyle),
                  dict(color='c',label="25mM NaCl",**CommonStyle),
                  dict(color='g',label="50mM NaCl",**CommonStyle),
                  dict(color='b',label="100mM NaCl",**CommonStyle)]
    # for style, just replace label 'NaCl' with 'KCL'
    StyleDictKCL = copy.deepcopy(StyleDicts)
    for i,TmpDict in enumerate(StyleDictKCL):
        TmpDict['label'] = TmpDict['label'].replace("NaCl","KCl")
        TmpDict['alpha'] = 0.7
    # determine the bounds for the FRET
    MinFret = -0.25
    MaxFretFromData = np.max([max(arr) for arr in data])
    MaxFret = 1.2
    # assume a-priori knowledge of the bin counts
    MinBin = 0
    MaxBin = 140
    StepFret = 0.01
    ylim = [MinBin,MaxBin]
    bins = np.arange(MinFret,MaxFret,StepFret)
    for i,histogram in enumerate(data):
        title = "Higher salt induces G-Quadruplex folding" \
                if i == 0 else ""
        # plot the NaCl data
        plt.subplot(NumHists,2,(2*i+1))
        plt.hist(histogram,bins=bins,**StyleDicts[i])
        AxesDefault(ylim,title,"Count",UseYTicks=True)
        plt.subplot(NumHists,2,2*(i+1))
        plt.hist(dataKcl[i],bins=bins,**StyleDictKCL[i])
        AxesDefault(ylim,"","",UseYTicks=False)
        # plot the KCl Data
    plt.subplot(NumHists,2,11)
    ax = plt.gca()
    plt.setp(ax.get_xticklabels(),visible=True)
    plt.setp(ax.get_yticklabels(),visible=True)
    pPlotUtil.lazyLabel("FRET","Count","",frameon=True)
    plt.subplot(NumHists,2,12)
    ax = plt.gca()
    plt.setp(ax.get_xticklabels(),visible=True)
    pPlotUtil.lazyLabel("FRET","","",frameon=True)

    pPlotUtil.savefig(fig,"./out.png")


if __name__ == "__main__":
    run()
