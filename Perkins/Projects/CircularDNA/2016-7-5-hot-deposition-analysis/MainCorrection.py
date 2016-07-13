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
    def __init__(self,MaxInvolsSizeMeters = 10e-9,FractionForOffset = 0.2,
                 SpatialGridUpSample = 5,MaxFourierSpaceComponent=10e-9):
        self.MaxInvolsSizeMeters = MaxInvolsSizeMeters
        self.FractionForOffset = FractionForOffset
        self.SpatialGridUpSample = SpatialGridUpSample
        self.MaxFourierSpaceComponent = MaxFourierSpaceComponent
    def CalculateOffset(self,Obj,Slice):
        """
        Calculates the force Offset for the given object and slice
        """
        return np.median(Obj.Force[Slice])
    def ZeroForceAndSeparation(self,Obj,IsApproach):
        """
        Given an object and whether it is the approach or retract, zeros it
        out 

        Args:
            Obj: TimeSepForce Object
            IsApproach: True if this represents the approach portion of the 
            Curve.
        Returns:
            Tuple of <Zeroed Separation, Zeroed Force>
        """
        Separation = Obj.Separation
        Time = Obj.Time
        SeparationZeroed = Separation - min(Separation)
        N = SeparationZeroed.size
        NOffsetPoints = int(np.ceil(N * self.FractionForOffset))
        # approach's zero is at the beginning (far from the surface)
        # retract is at the end (ibid)
        if (IsApproach):
            SliceV = slice(0,NOffsetPoints,1)
        else:
            SliceV = slice(-NOffsetPoints,None,1)
        # get the actual offset assocaiated with this object
        Offset = self.CalculateOffset(Obj,SliceV)
        return SeparationZeroed.copy(),(Obj.Force - Offset).copy()
    def FitInvols(self,Obj):
        """
        Fit to the invols on the (approach!) portion of Obj

        Args:
            Obj: TimeSepForceObject. We get just the approach from it and
            fit to that
        Returns:
            Nothing, but sets the object for future predicts
        """
        Approach,Retract = FEC_Util.GetApproachRetract(Obj)
        # get the zeroed force and separation
        SeparationZeroed,ForceZeroed  = self.ZeroForceAndSeparation(Approach,
                                                                    True)
        ArbOffset = max(np.abs(ForceZeroed))
        A = max(ForceZeroed)
        # adding in the arbitrary offset actually helps quite a bit.
        # we fit versus time, which also helps.
        FittingFunction = lambda t,tau :  np.log(A * np.exp(-t/tau)+ArbOffset)
        # for fitting, flip time around
        MaxTau = self.MaxInvolsSizeMeters
        params,_,_ = GenFit(SeparationZeroed,np.log(ForceZeroed+ArbOffset),
                            model=FittingFunction,
                            bounds=(0,MaxTau))
        # tau is the first (only) parameter
        self.Lambda= params[0]
        self.MaxForceForDecay = max(ForceZeroed)
    def PredictInvols(self,Obj,IsApproach):
        """
        Given an object, predicts the invols portion of its curve. *must* call
        after a fit

        Args:
            Obj: see FitInvols, except this is *either* the approach
            or retract
            IsApproach: see FitInvols
        Returns:
            Predicted, Zero-offset invols decay for Obj
        """
        SeparationZeroed,_ = self.ZeroForceAndSeparation(Obj,IsApproach)
        return self.MaxForceForDecay * np.exp(-SeparationZeroed/self.Lambda)
    def FitInterference(self,Obj):
        Approach,_ = FEC_Util.GetApproachRetract(Obj)
        # get the zeroed force and separation
        SeparationZeroed,ForceZeroed  = self.ZeroForceAndSeparation(Approach,
                                                                    True)
        # get the residuals (essentially, no 'invols') part
        FourierComponents = max(SeparationZeroed)/self.MaxFourierSpaceComponent
        NumFourierTerms = np.ceil(FourierComponents/self.SpatialGridUpSample)
        # down-spample the number of terms to match the grid
        # get the fourier transform in *space*. Need to interpolate onto
        # uniform gridding
        N = SeparationZeroed.size
        self.linear_grid = np.linspace(0,max(SeparationZeroed),
                                       N*self.SpatialGridUpSample)
        # how many actual terms does that translate into?
        ForceInterp =interp1d(x=SeparationZeroed,
                              y=Approach.Force,kind='linear')
        self.fft_coeffs = rfft(ForceInterp(self.linear_grid))
        # remove all the high-frequecy stuff
        NumTotalTermsPlusDC = int(2*NumFourierTerms+1)
        self.fft_coeffs[NumTotalTermsPlusDC:] = 0 
    def PredictInterference(self,Obj,IsApproach):
        # interpolate back to the original grid
        SeparationZeroed,_  = self.ZeroForceAndSeparation(Obj,
                                                          IsApproach)
        N = SeparationZeroed.size
        fft_representation = irfft(self.fft_coeffs)
        fft_pred = interp1d(x=self.linear_grid,
                            y=fft_representation)(SeparationZeroed)
        return fft_pred
    def CorrectApproachAndRetract(self,Obj):
        """
        Given an object, corrects and returns the approach and retract
        portions of the curve (dwell excepted)
        """
        Approach,Retract = FEC_Util.GetApproachRetract(Obj)
        SeparationZeroed,ForceZeroed = self.\
                    ZeroForceAndSeparation(Approach,IsApproach=True)
        # fit the invols
        self.FitInvols(Obj)
        # predict the invols -- we can subtract this off. 
        InvolsPrediction = self.PredictInvols(Approach,IsApproach=True)
        Approach.Force -= InvolsPrediction
        # now correct the interference artifct 
        self.FitInterference(Approach)
        fft_pred = self.PredictInterference(Approach,
                                            IsApproach=True)
        # make a prediction without the wiggles
        Approach.Force -= fft_pred
        # just for clarities sake, the approach has now been corrected
        ApproachCorrected = Approach
        InvolsPredictionRetract = self.PredictInvols(Retract,
                                                              IsApproach=False)
        RetractNoInvols = copy.deepcopy(Retract)
        RetractNoInvols.Force -= InvolsPredictionRetract
        # now correct the FFT stuff 
        fft_pred_retract = self.PredictInterference(RetractNoInvols,
                                                    IsApproach=False)
        RetractCorrected = copy.deepcopy(RetractNoInvols)
        RetractCorrected.Force -= fft_pred_retract
        return ApproachCorrected,RetractCorrected

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
    set_y_lim = lambda :  plt.ylim([-50,200])
    set_x_lim = lambda :  plt.xlim([-10,1300])
    PlotOptions = dict(PreProcess=True,
                       LegendOpts=dict(loc='upper left',frameon=True))
    for i,Tmp in enumerate(DataArray):
        Approach,Retract = FEC_Util.GetApproachRetract(Tmp)
        CorrectionObj = CorrectionObject()
        ApproachCorrected,RetractCorrected = \
                CorrectionObj.CorrectApproachAndRetract(Tmp)
        fig = pPlotUtil.figure()
        plt.subplot(2,1,1)
        FEC_Plot.FEC_AlreadySplit(Approach,Retract,XLabel="",**PlotOptions)
        set_y_lim()
        set_x_lim()
        plt.subplot(2,1,2)
        FEC_Plot.FEC_AlreadySplit(ApproachCorrected,RetractCorrected,
                                  **PlotOptions)
        plt.axhline(65,label="65pN",linewidth=3.0,linestyle='--')
        set_y_lim()
        set_x_lim()
        pPlotUtil.savefig(fig,"./tmp{:d}.png".format(i))

if __name__ == "__main__":
    run()
