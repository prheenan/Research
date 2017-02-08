# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy import signal
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities

def add_stretch(x,f,start,end,spring):
    f_copy = f.copy()
    idx_stretch = np.where( (x > start) & 
                            (x < end))
    x_stretch = x[idx_stretch]
    f_stretch = (x_stretch - x_stretch[0])**2 * spring
    return idx_stretch,f_stretch

def make_force_extension_curve(x,array_of_stretch_kwargs,DecayConst,snr):
    f = (1-np.exp(-x/DecayConst))
    for kwargs in array_of_stretch_kwargs:
        idx,f_at_idx = add_stretch(x=x,f=f,**kwargs)
        f[idx] += f_at_idx
    # add in uniform noise
    noise_ampltude = np.sqrt(1/snr)
    f += (np.random.normal(size=f.size)-0.5)* 2 * noise_ampltude
    return f

def normalized_force_extension(max_x=1,decay_const=1/50,points=1024,
                               spring_const=10,snr=100):
    x,dx = np.linspace(0,max_x,points,retstep=True)
    SpringStretch = spring_const
    DecayConst = decay_const
    stretch_kwargs = [
        # adhesion peak...
        dict(start=0.04,end=0.1,spring=SpringStretch*50),
        # next, add in something that looks (vaguely) like a WLC strech
        dict(start=0.2,end=0.4,spring=SpringStretch),
        dict(start=0.35,end=0.6,spring=SpringStretch*1.2),
        dict(start=0.55,end=0.9,spring=SpringStretch/2)]
    # get the force
    f = make_force_extension_curve(x,stretch_kwargs,DecayConst,snr)
    # convert to 0-1 (ish)
    f -= np.median(f)
    f /= np.max(f)
    # get the maximum force in 
    return x,dx,stretch_kwargs,f

def MakeFigure(Points=1024,MaxX=1,Seed=42,DecayConst=1/50,SpringStretch=10,
               snr=100):
    np.random.seed(Seed)
    x,dx,stretch_kwargs,f = \
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
