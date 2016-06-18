# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../../..")
from GeneralUtil.python import PlotUtilities as pPlotUtil



def run():
    """
    From Supplemental figure S1a 
    Vafabakhsh, R., and Ha, T. 2012) 
    Extreme Bendability of DNA Less than 100 Base Pairs Long Revealed by 
    Single-Molecule Cyclization. 
    """
    OverhangLength = [8,9,10]
    # we care about overhangs with 12...
    MaxLength = 12
    LoopingRatesPerMin = [0.1,0.05,0.08]
    UnloopingRatesPerMin = [0.4,0.04,0.008]
    # get a rough model for unlooping as a function of nucleotide
    coeffs = np.polyfit(x=OverhangLength,
                        y=np.log(UnloopingRatesPerMin),deg=1)
    xinterp = np.linspace(OverhangLength[0],MaxLength,num=100)
    Amp = coeffs[1]
    Decay = coeffs[0]
    PredFunc = lambda x : np.exp(Amp)*np.exp(x * Decay)
    PredRange = PredFunc(xinterp)
    # predict the looping rate (hopefully constant?)
    LoopingRateMean = np.mean(LoopingRatesPerMin)
    fig = pPlotUtil.figure(figsize=(8,12))
    plt.subplot(2,1,1)
    plt.semilogy(OverhangLength,LoopingRatesPerMin,'ko',label="Looping")
    plt.semilogy(OverhangLength,UnloopingRatesPerMin,'ws',label="Unlooping")
    plt.semilogy(xinterp,PredRange,'g-',label="Exponential fit to unlooping")
    plt.axhline(LoopingRateMean,linewidth=3,linestyle='--',
                label="Mean Looping Rate")
    pPlotUtil.lazyLabel("","Rate (min^(-1))",
                        "Vafabakhsh 2012: 12 nt overhang is >99.8% looped",
                        frameon=True)
    plt.xlim([7.5,12.5])
    plt.ylim([5e-5,1])
    plt.subplot(2,1,2)
    Looped = 1-PredRange/LoopingRateMean
    plt.plot(xinterp,Looped,
                 label="Predicted looped population")
    NumberLooped = Looped[-1]
    plt.plot(MaxLength,NumberLooped,'ro',
             label="{:.2f}% looped for 12nt".format(NumberLooped*100))
    plt.xlim([8.5,12.5])
    plt.ylim([0.5,1.05])
    pPlotUtil.lazyLabel("Overhang Length (nt)","Number looped at equilibrium",
                        "")
    pPlotUtil.savefig(fig,"Looping.png")

if __name__ == "__main__":
    run()
