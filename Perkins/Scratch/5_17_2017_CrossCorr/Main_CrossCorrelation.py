# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from numpy import convolve as convolve
from numpy import correlate as correlate



def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    TwoPi = 2 * np.pi
    # make the three regions in x
    n = 1000
    x1 = np.linspace(0,TwoPi,n)
    x2 = np.linspace(TwoPi,2*TwoPi,n)
    x3 = np.linspace(2*TwoPi,3*TwoPi,n)
    # make  the three regions in y
    y1 = -(2*(x1[0]-x1) + np.sin(x1-x1[0]))
    y1 -= y1[-1]
    y2 = np.zeros(x2.size)
    y3 = -3*(x3-x3[0]) + np.sin((x3-x3[0]))
    # concatenate them all
    x = np.concatenate((x1,x2,x3))
    y = np.concatenate((y1,y2,y3))
    # normalize the x and y
    norm = lambda x: (x-min(x))/(max(x)-min(x))
    x = norm(x)
    y = norm(y)
    # add in noise
    amplitude = 0.2
    AddNoise = lambda x : (x+(np.random.rand(x.size)-0.5)*2*amplitude)
    YNoiseAdded = AddNoise(y)
    q75, q25 = np.percentile(YNoiseAdded, [75 ,25])
    iqr = q75 - q25
    med = np.median(YNoiseAdded)
    StdNorm = lambda  x: (x-med)/iqr
    Normalize = lambda x : StdNorm(AddNoise(x))
    YNoise = StdNorm(YNoiseAdded)
    # get an offset slice
    SliceFract = 0.02
    NShift = int(YNoise.size*SliceFract)
    YShifted = Normalize(y[NShift:])
    XShifted = x[NShift:]
    DeltaX = np.median(np.diff(x))
    # note: I *assume* that YNoise is larger than YShifted in order to get
    # the true time shift.
    Convolved = correlate(YShifted,YNoise, mode='full')
    PointsConvolved = np.arange(0,Convolved.size,dtype=np.float64)
    MaxPoints = YNoise.size
    PointsConvolved = MaxPoints - PointsConvolved +1
    MaxConvolved = int(PointsConvolved[np.argmax(Convolved)])
    NFullPoints = YShifted.size+YNoise.size-1
    xlim = lambda : plt.xlim([0,max(x/DeltaX)])
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(x/DeltaX,StdNorm(y),'r-',linewidth=3.0,label="Noiseless")
    plt.plot(x/DeltaX,YNoise,'k-',label="Noisy, Unshifted",alpha=0.3)
    plt.plot(XShifted/DeltaX,YShifted,b'-',label="Noisy,Shifted",alpha=0.3)
    plt.legend(loc='lower center')
    xlim()
    plt.subplot(2,1,2)
    plt.plot(PointsConvolved,Convolved,'r--',label="Convolution")
    plt.axvline(MaxConvolved,
                label="Expected Shift {:d}".format(MaxConvolved))
    plt.axvline(NShift,linestyle='--',
                label="Actual Shift {:d}".format(NShift))
    plt.legend(loc='lower center')
    xlim()
    plt.show()

if __name__ == "__main__":
    run()
