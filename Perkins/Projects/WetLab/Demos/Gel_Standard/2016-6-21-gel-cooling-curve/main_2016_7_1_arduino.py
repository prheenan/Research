# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities as pPlotUtil

def GetTauAndPred(x,y):
    logy = np.log(y)
    coeffs = np.polyfit(x=x,y=logy,deg=1)
    tau = abs(1/coeffs[0])
    pred = np.exp(np.polyval(coeffs,x))
    print(tau)
    return tau,pred

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    File = "./ArduinoData/TimeInMilliSecondsTemperatureDegreesCelcius.csv"
    Arr = np.genfromtxt(File,delimiter=",")
    TimeSeconds = Arr[:,0]/1000
    TemperatureCelcius = Arr[:,1]
    # fridge is at 4C
    CoolAmbientTemperatureCelcius=4
    EndHeat = 600
    EndCool = 9800
    SliceHeat=slice(0,EndHeat)
    SliceCool=slice(EndHeat,EndCool)
    SliceFinal=slice(EndCool,None)
    SlicesAndStyles = [ [SliceHeat,dict(label="Initial Heating",color='r')],
                        [SliceCool,dict(label="4C Cooling",color='b')],
                        [SliceFinal,dict(label="Separation",color='r')]
                    ]
    #fit an exponential to the decay
    DecayY = TemperatureCelcius[SliceCool]
    DecayTime = TimeSeconds[SliceCool]
    Tau,Pred = GetTauAndPred(DecayTime-DecayTime[0],
                             DecayY-CoolAmbientTemperatureCelcius)
    TauHours = Tau/3600
    Pred += CoolAmbientTemperatureCelcius
    ToPlottingTime = lambda x : x/3600
    fig = pPlotUtil.figure()
    plt.plot(ToPlottingTime(TimeSeconds),TemperatureCelcius,'r.')
    plt.plot(ToPlottingTime(DecayTime),Pred,'b--',linewidth=4,
             label=(r"Temp(time)$\propto e^{-t/\tau}$" + 
                    r",$\tau$={:.1f} hours".format(TauHours)))
    pPlotUtil.lazyLabel("Time (Hours)","Temperature (Celcius)","")
    pPlotUtil.savefig(fig,"ArduinoOut.png")


if __name__ == "__main__":
    run()
