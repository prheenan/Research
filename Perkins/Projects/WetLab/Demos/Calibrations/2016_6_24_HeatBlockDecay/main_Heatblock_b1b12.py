# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../")
import GeneralUtil.python.PlotUtilities as pPlotUtil


def run():
    """
    Data taken on heatblock, SN 00151 on 2016-6-24,
    with inner six wells filled with water (90% filled) 
    """
    # note: ambient temperature is 22C
    AmbientTemp = 22
    # recording minutes, seconds, and temperature in centigrade, in that
    # order
    MinutesSecondsTemperature  =\
                                 [ [0,0,80.5],
                                   [2,0,78],
                                   [4,0,75],
                                   [7,15,70],
                                   [11,20,62],
                                   [16,10,59],
                                   [19,04,55.5],
                                   [24,00,51],
                                 ]
    # get just the times and temperatures
    Times = np.array([ t[0]*60 + t[1] for t in MinutesSecondsTemperature])
    Temperature = np.array([t[-1] for t in MinutesSecondsTemperature])
    # get the temperature, relative to ambient
    TemperatureRel = Temperature - AmbientTemp
    TimeRel = Times - Times[0]
    # get the log relative temperature
    LogTemp = np.log(TemperatureRel)
    # fit to an exponential decay
    coeffs = np.polyfit(x=TimeRel,y=LogTemp,deg=1)
    pred_log = np.polyval(coeffs,TimeRel)
    pred_temp = np.exp(pred_log) + AmbientTemp
    # get the tau (decay constant) in seconds
    tau= -1/coeffs[0]
    tau_minutes = tau/60
    # plot the exponential decay
    fig = pPlotUtil.figure()
    TimeMinutes = Times/60
    plt.errorbar(TimeMinutes,Temperature,yerr=0.5,fmt='bo-')
    plt.plot(TimeMinutes,pred_temp,'r',linewidth=3,
             label=r"Exponential Decay, $\tau$={:.1f} minutes".\
             format(tau_minutes))
    plt.axhline(AmbientTemp,label="Ambient Temperature",
                color="g",linestyle="--",linewidth=3)
    pPlotUtil.lazyLabel("Time (minutes)",r"Temperature ($^\circ$C)","")
    pPlotUtil.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
