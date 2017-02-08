# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy import signal
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Personal.EventDetection.Docs.ToyGraphs import SimulationUtil


def MakeFigure(Points=1024,MaxX=1,Seed=42,DecayConst=1/50,SpringStretch=10,
               snr=100):
    np.random.seed(Seed)
    x,dx,stretch_kwargs,f = SimulationUtil.\
        normalized_force_extension(max_x=MaxX,decay_const=DecayConst,
                                   points=Points,spring_const=SpringStretch,
                                   snr=snr)
    # first, add in approach/retract
    # make everything nice and zero-meaned, max of 1
    f -= np.mean(f)
    f /= max(f)
    # Get the FFT of our function
    fft_coeffs = np.fft.rfft(f)
    fft_freq = np.fft.rfftfreq(n=f.size,d=dx)
    # get the CWT of our function 
    CoeffMax = 50
    CoeffMin = 1
    NCoeffs = 50
    widths = np.linspace(CoeffMin,CoeffMax,NCoeffs)
    cwt_coeffs = signal.cwt(data=f,wavelet=signal.ricker,widths=widths)
    # see how well an out-of-the-box scipy system can do 
    peak_idx_scipy = signal.find_peaks_cwt(vector=f, widths=widths)
    # plot the various transforms..
    plt.subplot(3,1,1)
    plt.plot(x,f)
    PlotUtilities.lazyLabel("x","f(x)","Fourier and Wavelet Transform of f(x)")
    plt.subplot(3,1,2)
    plt.plot(fft_freq[1:]/max(fft_freq),fft_coeffs[1:])
    xlabel = r"$\frac{\mathrm{Frequency}}{\mathrm{Frequency}_{\mathrm{max}}}$"
    PlotUtilities.lazyLabel(xlabel,
                            "Positive FFT coefficients","")
    plt.subplot(3,1,3)
    # XXX make heat map like ch 6...
    plt.imshow(cwt_coeffs, extent=[0, max(x), min(widths), max(widths)], 
               cmap='PRGn', aspect='auto',
               vmax=abs(cwt_coeffs).max(), vmin=-abs(cwt_coeffs).max())
    title = "Positive wavelet coefficients for Laplacian of Gaussian"
    PlotUtilities.lazyLabel("x","Frequency",title)

def run():
    """
    simple script to compare wavelet transform and fft on something that
    looks like our data
    """             
    fig = PlotUtilities.figure(figsize=(10,16))
    MakeFigure()
    PlotUtilities.savefig(fig,"out.png")

if __name__ == "__main__":
    run()
