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
    # shift *both*, in differnet ways
    SliceNoise = slice(0,None,1)
    YNoise = AddNoise(y)[SliceNoise]
    YShifted = AddNoise(y[NShift:])
    XShifted = x[NShift:]
    XNoise = x[SliceNoise]
    return XShifted,XNoise,YShifted,YNoise,NShift

def TestSpeed(Sizes=np.logspace(1,7),SizeNaiveCutoff=5e4):
    """
    Speed tests everything

    Args:
       Sizes: to test the FFT-based cross correlation
       SizeNaiveCutoff: when the naive cross correlation is cutoff
    """
    from timeit import Timer
    times = []
    timesNaive = []
    sizesNaive = []
    for j,size in enumerate(Sizes):
        # get the data and shift it
        x,y = GetSampleFEC(size)
        _,_,YShifted,YNoise,_ = \
                AddNoiseAndShift(x,y,SliceFract=0.1)
        t = Timer(lambda: \
                  NormalizedCorrelation(YNoise,YShifted,
                                        CorrelationFunc=FFTCrossCorrelation))
        # get the expected convolution times
        BestTime = t.timeit(number=5)
        times.append(BestTime)
        # see if we should record naive data
        if (size < SizeNaiveCutoff):
            t = Timer(lambda: \
                NormalizedCorrelation(YNoise,YShifted,
                                      CorrelationFunc=NaiveCrossCorrelation))
            NaiveTime = t.timeit(number=5)
            timesNaive.append(NaiveTime)
            sizesNaive.append(size)
    # remained is just for plotting the results
    interpX = lambda x: np.linspace(min(x),max(x),
                                    len(x)*100)
    size_naive_interp = interpX(sizesNaive)
    # fit a quadratic to the naive
    coeffs = np.polyfit(x=sizesNaive,y=timesNaive,deg=2)
    quad = np.polyval(coeffs,size_naive_interp)
    # fit nlogn to the fft
    mFunc = lambda x,c: c*x*np.log(x)
    pOpt,_ = curve_fit(mFunc,xdata=Sizes,ydata=times,p0=1)
    # interpolate along sizes
    size_fft_interp = interpX(Sizes)
    nLogN = mFunc(size_fft_interp,*pOpt)
    fig = plt.figure()
    plt.loglog(Sizes,times,'ro',label="FFT Convolution")
    plt.loglog(sizesNaive,timesNaive,'bx',label="Naive Convolution")    
    plt.loglog(size_naive_interp,quad,'k--',label="O(n^2)")
    plt.loglog(size_fft_interp,nLogN,'r-',label="O(nlog(n))")
    plt.ylabel("Time (seconds)")
    plt.xlabel("Input Size (Number of Points)")
    plt.title("Naive convolution is O(n^2)")
    plt.legend(loc='upper left')
    MaxT= max(max(timesNaive),max(times))
    MinT = min(min(timesNaive),min(times))
    plt.ylim([min(times),MaxT])
    fig.savefig("./OutTimeTrials.png")

def TestCorrectness(ShiftPercentages=np.linspace(0,1,num=10),
                    rtol=1e-2,atol=1e-2,Sizes=None):
    """
    Tests the algorithm's correctness on a variety of shifts and data sizes

    Throws an error if it breaks.
    
    Args:
         ShiftPercentages: list of fractions [0,1] to shift the data
         rtol: relative tolerance allowed in the perceived fractional shift
         atol: absolute tolerance allowed in the perceived fractional shift
         (2e-2 means a different of 2% is ok)
      
         Sizes: Sizes of the data to use. If none, defaults to a 4 OOM range
    """
    N = ShiftPercentages.size
    # default sizes
    if (Sizes is None):
        Sizes = np.logspace(3,4,num=ShiftPercentages.size,base=10,
                            dtype=np.uint64)
    # loop through each percent and shift
    Errors = np.zeros((ShiftPercentages.size,Sizes.size))
    print("Doing Correction Tests...")
    for i,pct in enumerate(ShiftPercentages):
        for j,size in enumerate(Sizes):
            # get the data and shift it
            x,y = GetSampleFEC(size)
            XShifted,XNoise,YShifted,YNoise,NShift = \
                AddNoiseAndShift(x,y,SliceFract=pct)
            toTest = [ ([YNoise,YShifted],NShift),
                       ([YShifted,YNoise],NShift)]
            for args,Shift in toTest:
                # get the expected convolution
                PointsConvolved,Convolved = NormalizedCorrelation(*args)
                MaxConvolved = int(PointsConvolved[np.argmax(Convolved)])
                NPointsTotal = len(y)
                Errors[i][j] = abs(MaxConvolved-NShift)/NPointsTotal
                print("Checking if shift ({:d}) is close to predicted ({:d})".\
                      format(MaxConvolved,Shift))
                # if that fails, check relative
                np.testing.assert_allclose(NShift/NPointsTotal,
                                           MaxConvolved/NPointsTotal,
                                           atol=atol,rtol=rtol)
    fig = plt.figure()
    Styles = [dict(color='r',marker='o',linestyle='--'),
              dict(color='g',marker='s',linestyle='-'),
              dict(color='b',marker='v',linestyle='-.'),
              dict(color='k',marker='*',linestyle='-',linewidth=4.0)]
    for j,s in enumerate(Sizes):
        sty = Styles[j % len(Styles)]
        plt.plot(ShiftPercentages,Errors[:,j],label="Size={:d}".format(s),
                 **sty)
    plt.title("Error in shift determination <0.1%")
    plt.legend()
    plt.ylabel("Error (Percentage of Total number of points)")
    plt.xlabel("Amount Shifted")
    fig.savefig("./ErrorBySizeAndPct.png")
            
    
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
    XLimits = [min(PointsConvolved),max(PointsConvolved)]
    xlim = lambda : plt.xlim(XLimits)
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
    plt.legend(loc='upper left')
    xlim()
    plt.subplot(2,1,2)
    plt.plot(PointsConvolved,Convolved,'r-',markersize=2,label="Convolution")
    plt.axvline(MaxConvolved,
                label="Expect Shift: {:d} Pts".format(MaxConvolved))
    plt.axvline(NShift,linestyle='--',
                label="Actual Shift: {:d} Pts".format(NShift))
    plt.xlabel("Time (au)")
    plt.ylabel("Convolution (normalized)")
    plt.legend(loc='upper left')
    xlim()
    plt.show()
    fig.savefig(out)


    

def run(RunSpeedTests=False,seed=42):
    """
    Runs the tests according to the flags

    Args:
        Run<xx>: Runs <xx> if true
        seed: Psuedo-random number generator seed
    """
    np.random.seed(seed)
    PlotExampleCorrelation()
    TestCorrectness()
    print("Passed correctness tests.")
    if (RunSpeedTests):
        TestSpeed()

if __name__ == "__main__":
    run()
