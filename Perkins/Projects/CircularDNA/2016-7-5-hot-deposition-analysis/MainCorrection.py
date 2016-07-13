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
from FitUtil.FitUtils.Python.FitUtil import GenFit

class CorrectionObject:
    def __init__(self,DecayConstants,FourierCoefficients):
        self.DecayConstants = DecayConstants
        self.FourierCoefficients = FourierCoefficients
    def Correct(self,Separation,Force):
        pass

def ReadInData(FullName):
    mObjs = FEC_Util.ReadInData(FullName)
    return mObjs

    
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
    # figure out the first time we are close to zero, as an upper bound for
    # tau (spatial decay).
    #Generally, we shouldnt have more than 40nm (~160pN@4pN/nm)
    # for any reasonable invols, for DNA experiments
    MaxInvolsSizeMeters = 40e-9
    # fraction of start of approach we use for fitting
    FractionForOffset = 0.1
    # amount to up-sample for getting a uniform grid in space (For the Fourier
    # series
    SpatialGridUpSample = 5
    # maximum spatial resolution of Fourier series. This is 1/(2*f_nyquist),
    # where f_nyquist is the 'frequency' (inverse spatial component) associated
    # with the wiggle correction. The lower this is, the higher-order
    # the frequency correction is. Too high, and you will pick up noise. 
    MaxFourierSpaceComponent = 20e-9
    for i,Tmp in enumerate(DataArray):
        Approach,Retract = FEC_Util.GetApproachRetract(Tmp)
        Force = Approach.Force
        Time = Approach.Time
        Separation = Approach.Separation
        N = Force.size
        NumForOffset = int(N*FractionForOffset)
        # get the y offset, 'true zero'
        offset = np.median(Force[:NumForOffset])
        ForceZeroed = Force-offset
        ArbOffset = max(np.abs(ForceZeroed))
        A = max(ForceZeroed)
        SeparationZeroed = Separation - min(Separation)
        # adding in the arbitrary offset actually helps quite a bit.
        # we fit versus time, which also helps.
        FittingFunction = lambda t,tau :  np.log(A * np.exp(-t/tau)+ArbOffset)
        # for fitting, flip time around
        TimeFit = Time[::-1]
        MaxTau = MaxInvolsSizeMeters
        params,_,pred = GenFit(SeparationZeroed,np.log(ForceZeroed+ArbOffset),
                               model=FittingFunction,
                               bounds=(0,MaxTau))
        # tau is the first (only) parameter
        tau = params[0]
        MaxForce = max(ForceZeroed)
        Prediction = A*np.exp(-SeparationZeroed/tau)
        # get the residuals (essentially, no 'invols') part
        FourierComponents = max(SeparationZeroed)/MaxFourierSpaceComponent
        NumFourierTerms = int(np.ceil(FourierComponents/SpatialGridUpSample))
        # down-spample the number of terms to match the grid
        # get the fourier transform in *space*. Need to interpolate onto
        # uniform gridding
        linear_grid = np.linspace(0,max(SeparationZeroed),N*SpatialGridUpSample)
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
