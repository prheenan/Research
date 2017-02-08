# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy import signal
sys.path.append("../../../")
from GeneralUtil.python import PlotUtilities


def MakePlot():
    Points = 5000
    MaxX = np.pi * 2
    XValues,dx = np.linspace(start=0,stop=MaxX,num=Points,retstep=True)
    num_steps = 3
    step = int(Points/num_steps)
    frequency_low = 1
    frequency_high = 10
    freqs = np.linspace( frequency_low, frequency_high,num=num_steps )
    Idx = [ slice(i*step,(i+1)*step,1) for i in range(step)]
    y_values = np.zeros(XValues.size)
    for f,idx_range in zip(freqs,Idx):
        x_tmp = XValues[idx_range]
        y_values[idx_range] = np.sin(2*np.pi*f*(x_tmp-x_tmp[0]))
    # Get the FFT of our function
    fft_coeffs = np.fft.rfft(y_values)
    fft_freq = np.fft.rfftfreq(n=y_values.size,d=dx)
    # get the CWT of our function 
    CoeffMax = frequency_high*5
    CoeffMin = 1/(frequency_low/5)
    NCoeffs = 50
    widths = np.linspace(CoeffMin,CoeffMax,NCoeffs)
    wavelet_signal = lambda n_points,width: \
            signal.morlet(M=n_points,w=width,complete=True,s=1.0)
    cwt_coeffs = signal.cwt(data=y_values,wavelet=signal.ricker,widths=widths)
    fudge = 0.5
    x_lim = np.array([0,max(XValues)])
    x_fudge = np.array([-fudge,fudge])
    plt.subplot(3,1,1)
    plt.plot(XValues,y_values)
    plt.xlim(x_lim+x_fudge)
    PlotUtilities.lazyLabel("Time","","")
    plt.subplot(3,1,2)
    plt.plot(fft_freq,fft_coeffs)
    PlotUtilities.lazyLabel("Frequency","FFT Coefficients","")
    plt.subplot(3,1,3)
    plt.imshow(cwt_coeffs, extent=[x_lim[0], x_lim[1], 
                                   min(widths), max(widths)], 
               cmap='PRGn', aspect='auto',
               vmax=abs(cwt_coeffs).max(), vmin=-abs(cwt_coeffs).max())
    plt.xlim(x_lim+x_fudge)
    PlotUtilities.lazyLabel("Time","Frequency","Morlet Coefficient Map")


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    fig = PlotUtilities.figure(figsize=(10,16))
    MakePlot()
    PlotUtilities.savefig(fig,"out.png")

if __name__ == "__main__":
    run()
