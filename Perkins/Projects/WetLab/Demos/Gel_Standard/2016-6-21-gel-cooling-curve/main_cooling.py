# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities as pPlotUtil

class MinSec:
    def __init__(self,minutes,seconds=0):
        self.minutes = minutes
        self.seconds = seconds

class MeltData:
    def __init__(self,MinSec,Temperature):
        self.MinSec = MinSec
        self.Temperature=  Temperature
    @property
    def TimeInSeconds(self):
        return self.MinSec.minutes * 60 + self.MinSec.seconds

def run():
    """
    Demonstrates protocol for dna preparation
    """
    # using 1 tick/C thermeter
    TempError = 1
    # probably a 10s delay for recording
    TimeError = 10
    # bath temperature kept at about 1C
    MinTemp = 1
    # offset to start of 100V, in minutes
    Offset = 3*60+31+10
    # Temperature Data from 6/21/2016, consist of [time,temperature]
    Data = [
        # 200V
        [MinSec(0),58],
        # 0V
        [MinSec(10),64],
        [MinSec(12,10),62],
        [MinSec(15,55),59],
        [MinSec(17,20),57],
        [MinSec(19,35),55],
        [MinSec(21,56),53],
        [MinSec(25,30),50],
        [MinSec(31,46),46],
        [MinSec(42,50),40],
        [MinSec(56,27),34],
        [MinSec(69,00),30],
        [MinSec(88,53),24.5],
        [MinSec(97,13),22.5],
        [MinSec(109,05),20],
        [MinSec(122,0),18],
        [MinSec(135,40),16],
        [MinSec(160,45),13],
        [MinSec(170,10),10.5],
        [MinSec(212,50),7.5],
        # 100V
        [MinSec(Offset+7,54),7],
        [MinSec(Offset+21,0),8.5],
        [MinSec(Offset+25,10),9],
        [MinSec(Offset+42,05),9],
        [MinSec(Offset+55,30),9.5],
        [MinSec(Offset+69,07),10],
        [MinSec(Offset+89,15),10.5],
        [MinSec(Offset+111,50),11],
        [MinSec(Offset+118,50),11.5],
    ]
    # Convert the Data to something a bit easier to work with
    Objs = [MeltData(o,temp) for o,temp in Data]
    # get just the time and temperautres
    temp = np.array([o.Temperature for o in Objs])
    time = np.array([o.TimeInSeconds for o in Objs])
    # fit just the cooling portion to a an exponential decay model
    CoolingStartIdx = 1
    CoolingEndIdx = 19
    CoolSlice = slice(CoolingStartIdx,CoolingEndIdx,1)
    CoolingTime = time[CoolSlice]
    CoolingTemp = temp[CoolSlice]
    LogTemp = np.log(CoolingTemp)
    TimeStart = CoolingTime[0]
    TimeRel = CoolingTime - TimeStart
    coeffs = np.polyfit(y=LogTemp,x=TimeRel,deg=1)
    # now predict the temperature 
    HoursPredLog = np.log10(max(TimeRel))
    PredictX = np.logspace(-1,HoursPredLog,num=500)
    PredLogTemp = np.polyval(coeffs,x=PredictX)
    PredTemp = np.exp(PredLogTemp)
    # plot it.
    fig = pPlotUtil.figure()
    time_hours = time/3600
    plt.errorbar(x=time_hours,y=temp,yerr=TempError,
                 fmt='k*-',markersize=10,label="Measurements")
    # get the decay constant (in seconds)
    tau = 1/coeffs[0]
    # convert to a decay constant in sensible units, get the absolute value
    tau_hours = abs(tau/3600)
    pred_hours = (PredictX+TimeStart)/3600
    plt.plot(pred_hours,PredTemp,linewidth=3,linestyle="--",color='g',
             label=(r"Temp(time)$\propto e^{-t/\tau}$" + 
                    r",$\tau$={:.1f} hours".format(tau_hours)))
    # plot the various regions
    Spans = [ [[0,10/60],r"200V,(353$\rightarrow$395)mA"],
              [[10/60,max(pred_hours)],"0V, 0mA"],
              [[max(pred_hours),max(time_hours)*1.1],
               r"100V, (63$\rightarrow$67)mA"]]
    styles = [ dict(color='r'),
               dict(color='b'),
               dict(color='r')]
    for i,(region,label) in enumerate(Spans):
        color_idx = i % len(styles)
        plt.axvspan(*region,label=label,alpha=0.3,**(styles[color_idx]))
    # plot a line at the melting temperature
    MeltingTemp = 38
    plt.axhline(MeltingTemp,
                label=r"Circular DNA $T_m$={:d}$^\circ$C".\
                format(MeltingTemp),
                linewidth=4)
    pPlotUtil.lazyLabel("Time (hours)",r"Temperature ($^\circ$C)",
    "Three-stage purification denatures primers, discourages dimerization",
                        frameon=True)
    plt.ylim([-3,max(temp)*1.05])
    plt.xlim([-1/60,max(time_hours)*1.05])
    pPlotUtil.savefig(fig,"./MT.png")

if __name__ == "__main__":
    run()
