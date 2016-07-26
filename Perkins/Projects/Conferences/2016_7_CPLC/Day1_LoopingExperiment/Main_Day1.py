# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../../")
from GeneralUtil.python import PlotUtilities as pPlotUtil
from FitUtil.FitUtils.Python import FitUtil as FitUtil


def LoopingModel(x,tau,a):
    """
    model for exponential increase to a value

    Args:
        x: indpendent parameter, same units as tau
        tau: decay constant
        a: scaling constant, grow to this
    Returns:
        see below
    """
    return a * ( 1-np.exp(-x/tau))
    
def BoundedFit(x,y,TauBounds):
    """
    Fit LoopingModel to x and y, given bounds on tau. The bounds on the 
    maximum are assumed to coem from y (cant go above the max)

    Args:
       x : see LoopingModel
       y : values we are tyring to fit to loopingModel
       TauBounds: bounds on the time constant
    Returns:
       tuple of <params,paramstdev>, where element 0 is tau, element 1
       is the amplitude
    """
    # amplitude can't be larger that y, should be positive
    AmplitudeBounds = [0,max(y)]
    # combine the bounds
    BoundsLow = [i[0] for i in [TauBounds,AmplitudeBounds]]
    BoundsHigh = [i[1] for i in [TauBounds,AmplitudeBounds]]
    Bounds = (BoundsLow,BoundsHigh)
    # write down the model to use
    # fit everything
    params,stdevs,pred = FitUtil.GenFit(np.array(x),np.array(y),
                                        LoopingModel,bounds=Bounds,
                                        # for initial guesses, just guess the
                                        # mean
                                        p0 = [np.mean(TauBounds),
                                              np.mean(AmplitudeBounds)])
    return params,stdevs

def GetPredictedXandY(x,params):
    """
    Given x and params, gets  the predicted x and y by LoopingModel

    Args:
        x: See LoopingModel
        y: params: first tuple of return of BoundedFit
    Returns:
        tuple of (inteprolated, higher res) x and params
    """
    InterpolatedX = np.linspace(min(x),max(x),np.array(x).size*10)
    PredictedY = LoopingModel(InterpolatedX,*params)
    return InterpolatedX,PredictedY

def GetExponentialFit(MinTau,MaxTau,x,y):
    """
    given x,y, and bounds for tau, Fits looping model, getting the fitted
    x and y, and associated params and stdevs
    
    Args:
       Min/MaxTau: lower an upper bounds for decay time cons
       x,y: see BoundedFit
    Returns
        tuple of Interpolated (10 spacing) X,Y,Params,ParamsStde
    """
    TauBounds = [MinTau,MaxTau]
    Params,ParamsStd = BoundedFit(x,y,TauBounds)
    FitX,FitY = GetPredictedXandY(x,Params)
    return FitX,FitY,Params,ParamsStd

def run():
    """
    Creates a graph of looping rate versus time for our fret data
    """
    # write down the numerator (number in high fret) and denominator
    # (total number) for the third sample, which was digested after
    # 1 minute in looping buffer
    Numer3 = [19,185,282,261,340,315,301,276,290,245]
    Denom3 = [642,473,538,428,450,424,381,365,371,303]
    off = 3
    Times3 = [0,off+0,off+4,off+8,off+17,off+23,off+31,off+35,off+44,off+53]
    # ibid, example fourth sample, which was 90 minutes
    Numer4 = [21,175,349,435,470,495,552,536,504,524,455]
    Denom4 = [922,811,906,869,809,728,757,696,663,674,639]
    offFour = off
    Times4 = [0,offFour+0,offFour+3,offFour+8,offFour+11,offFour+20,
              offFour+25,offFour+33,offFour+37,offFour+47,offFour+59]
    # get the population ratio
    Ratio3 = np.array(Numer3)/np.array(Denom3)
    Ratio4 = np.array(Numer4)/np.array(Denom4)
    # fit the exponential 'decays' 
    MaxTau = 10
    MinTau = 1
    ThreeX,ThreeY,ThreeParams,ThreeStd = GetExponentialFit(MinTau,MaxTau,
                                                           Times3,Ratio3)
    FourX,FourY,FourParams,FourStd = GetExponentialFit(MinTau,MaxTau,
                                                       Times4,Ratio4)
    # make our figure
    fig = pPlotUtil.figure()
    LabelThree = "1 minute cutoff library\n" + \
                 r"$\tau$={:.1f} minutes".format(ThreeParams[0])
    LabelFour =  "90 minute cutoff library\n" + \
                 r"$\tau$={:.1f} minutes".format(FourParams[0])
    plt.plot(Times3,Ratio3,'ro',label=LabelThree)
    plt.plot(Times4,Ratio4,'bs',label=LabelFour)
    plt.plot(ThreeX,ThreeY,'r-',linewidth=2)
    plt.plot(FourX,FourY,'b--',linewidth=2)
    MaxT = max(np.max(Times3),np.max(Times4))*1.01
    plt.xlim([-1,MaxT])
    pPlotUtil.lazyLabel("Time (minutes)","Population fraction in high fret",
                        "Lower cutoff library loops faster",frameon=True)

    pPlotUtil.savefig(fig,"Out.png")
    
if __name__ == "__main__":
    run()
