# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../..")
import GeneralUtil.python.PlotUtilities as pPlotUtil

class Cantilever:
    def __init__(self,Name,Length,Radius,Height):
        """
        Drag from 
        Mechanics of motor proteins and the cytoskeleton

        Args:
            Name,Length,Thickness,Height: of the cantilever (height relative
            to the surface
        """
        self.Name = Name
        self.L0 = Length
        self.h = Height
        self.r = Radius
    @property
    def FrictionPerEta(self):
        """
        Drag coefficient, per eta, from 
        Mechanics of motor proteins and the cytoskeleton
        """
        return 2*np.pi*self.L0/(np.arccosh(self.h/self.r))


def run():
    """
    Biolever long from http://probe.olympus-global.com/en/product/bl_rc150vb_hw/
    """
    Levers = [
        Cantilever(Name="Biolever A ('Short')",
                   Length=60e-6,
                   # radius as the thickness
                   Radius=0.18e-6,
                   # height from tip height
                   Height=7e-6),
        Cantilever(Name="Biolever B ('Long')",
                   Length=100e-6,
                   # radius as the thickness
                   Radius=0.18e-6,
                   # height from tip height
                   Height=7e-6),
        Cantilever(Name="Biolever Mini",
                   Length=38e-6,
                   # radius as half the width
                   Radius=0.2e-6,
                   # height from tip height
                   Height=7e-6),
        # See:
        # omega.albany.edu:8008/calc3/vector-functions-dir/dna-solution-m2h.html
        Cantilever(Name="DNA",
                   Length=650e-9,
                   Radius=1e-9,
                   Height=2e-9)
        ]
    # get the force in SI units
    Velocities = np.logspace(-9,-6)
    Forces = []
    for o in Levers:
        Forces.append(o.FrictionPerEta * Velocities)
    # Plot everything in pN and nm/s
    VelNmPerSec = Velocities * 1e9
    fig = pPlotUtil.figure(figsize=(8,8))
    ax = plt.subplot(1,1,1)
    for i,f in enumerate(Forces):
        ForcePn = f *1e12
        plt.plot(VelNmPerSec,ForcePn,label="Drag on " + Levers[i].Name)
    # write down the rupture force in pN based on
    """
    Hatch, K., Danilowicz, C., Coljee, V., and Prentiss, M. (2008).
    Demonstration that the shear force required to separate short 
    double-stranded DNA does not increase significantly with sequence length 
    for sequences longer than 25 base pairs. Phys. Rev. E 78, 011920.
    """
    RuptureForce = 20
    plt.axhline(RuptureForce,linewidth=3.0,linestyle="--",
                label="{:d}pN (Circular Rupture)".format(RuptureForce))
    TargetVel = 150
    plt.axvline(TargetVel,linewidth=3.0,linestyle="--",
                label="V={:d}nm/s".format(TargetVel))
    pPlotUtil.lazyLabel("Velocity (nm/s)","Force (pN)",
                        "V~1um/s in water likely to " +
                        "rupture circular DNA",loc='lower right',frameon=True)
    ax.set_yscale('log')
    pPlotUtil.savefig(fig,"DragForceAtVelocities.png")
    

if __name__ == "__main__":
    run()
