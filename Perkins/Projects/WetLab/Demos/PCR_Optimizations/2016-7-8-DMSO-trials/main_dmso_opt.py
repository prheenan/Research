# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities as pPlotUtil

class DmsoInfo:
    def __init__(self,DMSO,Repeats,YieldConcentration,YieldVolumeUl=35):
        """
        <Description>
        
        Args:
            param1: This is the first param.
    
        Returns:
            This is a description of what is returned.
        """
        self.DmsoPercent = DMSO
        self.RepeatNumber = Repeats
        self.YieldConc = YieldConcentration
        self.YieldVolume = YieldVolumeUl
        

def run():
    """
    """

    DMSO = [
            # data from 7/8/2016, all overhang spacers
            DmsoInfo(0,40,67),
            DmsoInfo(1,40,57),
            DmsoInfo(3,40,139),
            DmsoInfo(5,40,97)
    ]
    Pct = [d.DmsoPercent for d in DMSO]
    YieldsNanograms = [d.YieldConc*d.YieldVolume for d in DMSO]
    YieldsMicrograms = np.array(YieldsNanograms)/1000
    fig = pPlotUtil.figure()
    plt.plot(Pct,YieldsMicrograms,'ro')
    pPlotUtil.lazyLabel("DMSO %","Yield (ug)",
                        "Increasing DMSO Percentage increases yield")
    fudge = 0.5
    plt.xlim([-fudge,max(Pct)+fudge])
    pPlotUtil.savefig(fig,"DMSO.png")


if __name__ == "__main__":
    run()
