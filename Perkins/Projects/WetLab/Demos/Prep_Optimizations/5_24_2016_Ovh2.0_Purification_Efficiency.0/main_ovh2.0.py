# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

class Constructs:
    CIRCULAR_1p8KB_OVH2p0 = 0
    LINEAR_1p8KB = 1

class Method:
    FREEZE_SQUEEZE_OVERNIGHT = 0
    FREEZE_SQUEEZE_5MIN = 1

class Prep:
    def __init__(self,Construct,ConcLoaded,VolumeLoaded,AgarosePct,
                 YieldConc,YieldVol,
                 MethodV,Descr):
        """
        Args:
            Construct: class variable of Constructs
            ConcLoaded: Concentration loaded, in ng/uL
            VolumeLoaded: amount of concentraton loaded, in uL
            YieldConc: amount yielded, in ng/uL
            YieldVol: volume yielded, in ng/uL
            MethodV: the method used for purification
            Descr: string description
        """
        self.Construct = np.array(Construct)
        self.Concentration = ConcLoaded
        self.Volume = np.array(VolumeLoaded)
        self.AgarosePct = AgarosePct
        self.Method = Method
        self.MassLoadedNg = self.Concentration * self.Volume
        self.YieldConc = np.array(YieldConc)
        self.YieldVol = np.array(YieldVol)
        self.MassObtained = (self.YieldConc * self.YieldVol)
        self.Decr = Descr
    @property
    def Efficiency(self):
        """
        Returns actual mass obtained over loaded mass
        """
        return np.array(self.MassObtained)/np.array(self.MassLoadedNg)

def run():
    """
    Creates plots describing how the ovh2.0-idspacer works 
    """
    mObj = [ \
             # 6/6/2016
             Prep(Construct=Constructs.CIRCULAR_1p8KB_OVH2p0,
                  # note: all volumes in uL, all masses in ng
                  ConcLoaded=50,
                  VolumeLoaded=190,
                  # XXX fix...
                  YieldConc=[14,13,11],
                  AgarosePct=2,
                  YieldVol=[35,35,35],
                  MethodV=Method.FREEZE_SQUEEZE_OVERNIGHT,
                  Descr="Ovh2.0 Circular Band"),
             # 6/7/2016
             Prep(Construct=Constructs.CIRCULAR_1p8KB_OVH2p0,
                  # note: all volumes in uL, all masses in ng
                  ConcLoaded=44,
                  # loaded 2*190 gels
                  VolumeLoaded=190*2,
                  # XXX fix...
                  YieldConc=[31.5],
                  AgarosePct=1,
                  YieldVol=[42],
                  MethodV=Method.FREEZE_SQUEEZE_OVERNIGHT,
                  Descr="Ovh2.0 Circular Band"),
             # 6/9/2016
             Prep(Construct=Constructs.CIRCULAR_1p8KB_OVH2p0,
                  # note: all volumes in uL, all masses in ng
                  ConcLoaded=46,
                  # loaded 2*190 gels
                  VolumeLoaded=190,
                  # XXX fix...
                  YieldConc=[54,61.6],
                  AgarosePct=1,
                  YieldVol=[35],
                  MethodV=Method.FREEZE_SQUEEZE_5MIN,
                  Descr="Ovh2.0 Circular Band"),
             ]
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    for i,o in enumerate(mObj):
        Size = np.ones(o.Efficiency.size) * o.AgarosePct
        label = "Data" if i == 0 else ""
        plt.plot(Size,o.Efficiency,'ro',label=label)
    plt.axhline(0.2,linestyle='--',color='r',
                label="0.2 efficiency would be awesome")
    plt.xlabel("Agarose %")
    plt.ylabel("Efficiency [0,1]")
    plt.title("Efficiency of DNA prepraration as a function of Agarose %")
    plt.legend()
    plt.ylim([0,0.4])
    plt.xlim([0,2.5])
    fig.savefig("Labelled DNA efficiency")

if __name__ == "__main__":
    run()
