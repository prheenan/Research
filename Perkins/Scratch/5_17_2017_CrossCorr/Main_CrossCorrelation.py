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
    PointsConvolved = MaxPoints - PointsConvolved +1
    return PointsConvolved,Convolved

def PlotExampleCorrelation(n=50000,out="./ExampleCorr.png"):
    """
    Plots an example correlation

    Args:
        n: number of points to use
        out: where to save the result
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
    # add in noise
    amplitude = 0.1
    Range = max(y)-min(y)
    AddNoise = lambda x : (x+(np.random.rand(x.size)-0.5)*2*amplitude*Range)
    YNoise = AddNoise(y)
    # get an offset slice
    SliceFract = 0.02
    NShift = int(YNoise.size*SliceFract)
    YShifted = AddNoise(y[NShift:])
    XShifted = x[NShift:]
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
    plt.plot(x/DeltaX,StdNorm(YNoise,*NormBy),
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
    PlotExampleCorrelation()

if __name__ == "__main__":
    run()
