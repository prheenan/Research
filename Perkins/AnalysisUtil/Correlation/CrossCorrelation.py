# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from numpy import convolve as convolve
from numpy import correlate as correlate
from scipy.signal import fftconvolve
# curve fitting, for nlogn
from scipy.optimize import curve_fit


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

def NaiveCrossCorrelation(y1,y2):
    """
    Uses the (non FFT, O(n^2)) naive cross correlation to get the convolution

    Args:
         y1,y2 : find the convolution from y1 to y2
    Returns:
         normalized cross-correlation
    """
    Convolved = correlate(y1,y2, mode='full')
    return Convolved

def FFTCrossCorrelation(y1,y2):
    """
    Uses FFT convolution for an O(n*log(n)) algorithm for cross correlation 

    Args:
        y1,y2: see NaiveCrossCorrelation
    Returns:
        normalized cross-correlation
    """
    # correlation is just a convolution of the two arrays, one reversed.
    # see : https://en.wikipedia.org/wiki/Cross-correlation
    # also:
# stackoverflow.com/questions/12323959/fast-cross-correlation-method-in-python
    return fftconvolve(y1, y2[::-1], mode='full')
    

def NormalizedCorrelation(y1,y2,CorrelationFunc=FFTCrossCorrelation):
    """
    Returns offset information from the larger of  (y1,y2) to the smaller 
    of (y1,y2), defaulting to y1 to y2 if they tie. return is in units of data 
    points, according to the highest cross correlation.

    Args:
        y1,y2: the two arrays to check.
    Returns:
        Tuple of <list of indices from y1 to y2, normalized convolution 
        amplitudes associated>
    """
    # normalize the x and y
    norm = lambda x: (x-min(x))/(max(x)-min(x))
    # may need to swap y1 and y2 (we assume y1 is the larger)
    if (y2.size > y1.size):
        tmp = y1.copy()
        y1 = y2
        y2 = tmp
    YNoiseAdded = y1
    Normalize = lambda x : StdNorm(x,med,iqr)
    YNoise = StdNorm(y1,*Stats(y1))
    YShifted = StdNorm(y2,*Stats(y1))
    # note: I *assume* that YNoise is larger than YShifted in order to get
    # the true time shift.
    Convolved = CorrelationFunc(YShifted,YNoise)
    # normalize the convolved function
    MinV = min(Convolved)
    MaxV= max(Convolved)
    Convolved = (Convolved - MinV)/(MaxV-MinV)
    # determine the actual 'shift' necessary
    # XXX explicitly define this
    PointsConvolved = np.arange(0,Convolved.size,dtype=np.float64)
    MaxPoints = max(YNoise.size,YShifted.size)
    PointsConvolved = (MaxPoints - PointsConvolved -1)
    return PointsConvolved,Convolved
