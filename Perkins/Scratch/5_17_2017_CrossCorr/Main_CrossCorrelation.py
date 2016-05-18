# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from numpy import convolve as convolve
from numpy import correlate as correlate

def StdNorm(ToNorm,Offset,Scale):
    """
    Normalize an array to an offset and scale

    Args:
        ToNorm: What to normalize
        Offset: what the new zero is. commonly mean
        Scale: what the new scaling is. Commonly a standadrd deviation
    Returns: 
         The scaled array
    """
    return (ToNorm-Offset)/Scale

def Stats(x):
    """
    Returns:
        the median and IQR (75%-25%)
    """
    q75, q25 = np.percentile(x, [75 ,25])
    iqr = q75 - q25
    med = np.median(x)
    return med,iqr

def NormalizedCorrelation(y1,y2):
    """
    Returns the offset *from* y1 to y2, in units of data points, according to
    the highest cross correlation.

    Args:
        y1,y2: the two arrays to check.
    Returns:
        Tuple of <list of indices from y1 to y2, normalized convolution 
        amplitudes associated>
    """
    # normalize the x and y
    norm = lambda x: (x-min(x))/(max(x)-min(x))
    YNoiseAdded = y1
    Normalize = lambda x : StdNorm(x,med,iqr)
    YNoise = StdNorm(y1,*Stats(y1))
    YShifted = StdNorm(y2,*Stats(y1))
    # note: I *assume* that YNoise is larger than YShifted in order to get
    # the true time shift.
    Convolved = correlate(YShifted,YNoise, mode='full')
    Convolved = (Convolved - min(Convolved))/(max(Convolved)-min(Convolved))
    PointsConvolved = np.arange(0,Convolved.size,dtype=np.float64)
    MaxPoints = YNoise.size
    PointsConvolved = MaxPoints - PointsConvolved -1
    return PointsConvolved,Convolved

def GetSampleFEC(n):
    """
    Get an AFM-style Force Extensioon Curve (FEC)

    Args:
         n: number of points
    Returns:
         tuple of <x,y>, where x is time and y is the 'force'
    """
    TwoPi = 2 * np.pi
    # make the three regions in x
    x1 = np.linspace(0,TwoPi,n)
    x2 = np.linspace(TwoPi,2*TwoPi,n)
    x3 = np.linspace(2*TwoPi,3*TwoPi,n)
    ## make  the three regions in y
    # y1: this should 'grow' to zero ('approach')
    y1 = -(2*(x1[0]-x1) + np.sin(x1-x1[0]))
    y1 -= y1[-1]
    y2 = np.zeros(x2.size)
    # use the offset x to get something that 'decays' (invols-like) from 0s
    y3 = -3*(x3-x3[0]) + np.sin((5*x3-x3[0]))  + np.sin((50*x3-x3[0])) 
    # concatenate them all
    x = np.concatenate((x1,x2,x3))
    y = np.concatenate((y1,y2,y3))
    return x,y

def AddNoiseAndShift(x,y,SliceFract=0.02,amplitude=0.01):
    """
    Given x and y, adds amplitude 'amplitude', which is a fraction
    (ie: 0. to 1.) of max(y)-min(y)

    Args:
        x: x values, will be shifted
        y: y values, will have noise added and shifted
        SliceFract: how much to shift y by, in [0,1] fraction of len(y)
        amplitude: noise to add,  is a fraction 
        (ie: 0. to 1.) of max(y)-min(y)
    Returns:
        tuple of <XShifted,XNoise,YShifted,YNoise,NShift>
        where N shift is the actual number of points shifted
    """
    # add in noise
    amplitude = 0.1
    Range = max(y)-min(y)
    scale = amplitude * Range
    AddNoise = lambda x : (x+np.random.normal(loc=0,
                                              scale=scale,
                                              size=x.size))
    # get an offset slice
    SliceFract = 0.02
    NShift = int(y.size*SliceFract)
    YNoise = AddNoise(y)[:]
    YShifted = AddNoise(y[NShift:])
    XShifted = x[NShift:]
    XNoise = x[:]
    return XShifted,XNoise,YShifted,YNoise,NShift


def TestCorrectness(ShiftPercentages=np.linspace(0,1),rtol=1e-2,atol=1e-2,
                    AbsoluteShiftAllowed=4,Sizes=None):
    """
    Tests the algorithm's correctness on a variety of shifts and data sizes

    Throws an error if it breaks.
    
    Args:
         ShiftPercentages: list of fractions [0,1] to shift the data
         rtol: relative tolerance allowed in the perceived fractional shift
         atol: absolute tolerance allowed in the perceived fractional shift
         (2e-2 means a different of 2% is ok)
      
         Sizes: Sizes of the data to use. If none, defaults to a 4 OOM range
         AbsoluteShiftAllowed: how many data points is OK, from an absolute
         point of view
    """
    N = ShiftPercentages.size
    # default sizes
    if (Sizes is None):
        Sizes = np.logspace(1,4,num=ShiftPercentages.size,base=10)
    # loop through each percent and shift
    for pct in ShiftPercentages:
        for size in Sizes:
            # get the data and shift it
            x,y = GetSampleFEC(size)
            XShifted,XNoise,YShifted,YNoise,NShift = \
                AddNoiseAndShift(x,y,SliceFract=pct)
            # get the expected convolution
            PointsConvolved,Convolved = NormalizedCorrelation(YNoise,YShifted)
            MaxConvolved = int(PointsConvolved[np.argmax(Convolved)])
            print(MaxConvolved,NShift,len(y))
            NPointsTotal = len(y)
            # first,  check absolute
            if (np.allclose(NShift,MaxConvolved,atol=AbsoluteShiftAllowed,
                            rtol=0)):
                continue
            # if that fails, check relative
            np.testing.assert_allclose(NShift/NPointsTotal,
                                       MaxConvolved/NPointsTotal,
                                       atol=atol,rtol=rtol)

            
    
def PlotExampleCorrelation(n=10000,out="./ExampleCorr.png",**kwargs):
    """
    Plots an example correlation, saving to 'out'

    Args:
        n: number of points to use
        out: where to save the result
        kwargs: passed to AddNoiseAndShift after x and y
    """
    x,y = GetSampleFEC(n)
    XShifted,XNoise,YShifted,YNoise,NShift = AddNoiseAndShift(x,y,**kwargs)
    PointsConvolved,Convolved = NormalizedCorrelation(YNoise,YShifted)
    MaxConvolved = int(PointsConvolved[np.argmax(Convolved)])
    NFullPoints = YShifted.size+YNoise.size-1
    DeltaX = np.median(np.diff(x))
    xlim = lambda : plt.xlim([0,max(x/DeltaX)])
    # normalize eveythning
    NormBy = Stats(y)
    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.title("Correlation by Convolution efficiently determines time shift")
    plt.plot(x/DeltaX,StdNorm(y,*NormBy),'r-',linewidth=3.0,label="Noiseless")
    plt.plot(XNoise/DeltaX,StdNorm(YNoise,*NormBy),
             'k-',label="Noisy, Unshifted",alpha=0.3)
    plt.plot(XShifted/DeltaX,StdNorm(YShifted,*NormBy),b'-',
             label="Noisy,Shifted",alpha=0.3)
    plt.ylabel("Normalized measurement.")
    plt.legend(loc='lower center')
    xlim()
    plt.subplot(2,1,2)
    plt.plot(PointsConvolved,Convolved,'r--',label="Convolution")
    plt.axvline(MaxConvolved,
                label="Expect Shift: {:d} Points".format(MaxConvolved))
    plt.axvline(NShift,linestyle='--',
                label="Actual Shift: {:d} Points".format(NShift))
    plt.xlabel("Time (au)")
    plt.ylabel("Convolution (normalized)")
    plt.legend(loc='upper center')
    xlim()
    fig.savefig(out)


    

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    np.random.seed(42)
    PlotExampleCorrelation()
    TestCorrectness()
        

if __name__ == "__main__":
    run()
