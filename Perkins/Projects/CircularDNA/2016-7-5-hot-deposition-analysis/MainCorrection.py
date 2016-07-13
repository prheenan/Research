# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis.DataCorrection.\
    CorrectionMethods import CorrectForcePullByMetaInformation

from FitUtil.WormLikeChain.Python.Code.WLC_Fit import BoundedWlcFit
from FitUtil.FitUtils.Python.FitClasses import GetBoundsDict
from GeneralUtil.python import PlotUtilities as pPlotUtil
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
import copy
from GeneralUtil.python.IgorUtil import SavitskyFilter
from scipy.fftpack import rfft,irfft
from scipy.interpolate import interp1d

class CorrectionObject:
    def __init__(self,DecayConstants,FourierCoefficients):
        self.DecayConstants = DecayConstants
        self.FourierCoefficients = FourierCoefficients
    def Correct(self,Separation,Force):
        pass

def ReadInData(FullName):
    mObjs = FEC_Util.ReadInData(FullName)
    return mObjs

def GetFittingObject(Obj):

    return CorrectionObject()
    
def run():
    """
    Runs contour length analysis
    """
    OutFile = ""
    Limit = 2
    FullNames = ["2016_7_10_1ng_ul_50C_4hour_depo_circ_dna_Strept_tip_I.pxp"]
    DataArray = []
    for i,Name in enumerate(FullNames):
        DataArray.extend(pCheckUtil.getCheckpoint("Tmp{:d}.pkl".format(i),
                                                  ReadInData,False,Name))
    for i,Tmp in enumerate(DataArray):
        Approach,Retract = FEC_Util.GetApproachRetract(Tmp)
        Force = Approach.Force
        Separation = Approach.Separation
        N = Force.size
        NumForOffset = int(N/10)
        offset = np.median(Force[:NumForOffset])
        ForceZeroed = Force-offset
        SeparationZeroed = Separation - min(Separation)
        MaxForce = max(ForceZeroed)
        # XXX hand-picked exponential for now...
        Prediction = MaxForce*np.exp(-SeparationZeroed/1e-9)
        # get the residuals (essentially, no 'invols') part
        MaxFourierSpaceComponent = 100e-9
        NumFourierTerms = int(max(SeparationZeroed)/MaxFourierSpaceComponent)
        # get the fourier transform in *space*. Need to interpolate onto
        # uniform gridding
        linear_grid = np.linspace(0,max(SeparationZeroed),N*10)
        # how many actual terms does that translate into?
        ForceWithoutInvols = ForceZeroed-Prediction
        ForceInterp =interp1d(x=SeparationZeroed,
                              y=ForceWithoutInvols,kind='linear')
        fft_coeffs = rfft(ForceInterp(linear_grid))
        # remove all the high-frequecy stuff
        fft_coeffs[2*NumFourierTerms+1:] = 0 
        fft_representation = irfft(fft_coeffs)
        # interpolate back to the original grid
        fft_pred = interp1d(x=linear_grid,
                            y=fft_representation)(SeparationZeroed)
        # make a prediction without the wiggles
        WithoutWiggles = ForceWithoutInvols - fft_pred
        fig = pPlotUtil.figure(figsize=(6,6))
        ToXPlot = lambda x: x * 1e6
        ToYPlot = lambda y: y * 1e12
        NumPlots = 3
        y_lim_fixed = [-20,20]
        plt.subplot(NumPlots,1,1)
        plt.plot(ToXPlot(SeparationZeroed),ToYPlot(ForceZeroed),
                 color='k',alpha=0.3)
        plt.plot(ToXPlot(SeparationZeroed),ToYPlot(Prediction),
                 linewidth=4,linestyle='--')
        pPlotUtil.lazyLabel("","Force","")
        plt.subplot(NumPlots,1,2)
        plt.plot(ToXPlot(SeparationZeroed),ToYPlot(ForceWithoutInvols))
        plt.plot(ToXPlot(SeparationZeroed),ToYPlot(fft_pred))
        plt.ylim(y_lim_fixed)
        pPlotUtil.lazyLabel("","Force","")
        plt.subplot(NumPlots,1,3)
        plt.plot(ToXPlot(SeparationZeroed),ToYPlot(WithoutWiggles))
        pPlotUtil.lazyLabel("Separation","Force","")
        plt.ylim(y_lim_fixed)
        pPlotUtil.savefig(fig,"out{:d}.png".format(i))
        

if __name__ == "__main__":
    run()
