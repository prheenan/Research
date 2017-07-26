# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


sys.path.append("../../../../")
from Util import DilutionUtil 

def run():
    print("No Divalent buffer creation, 2x")
    Stats = [ ["Tris-HCl","mM",1000,20,0],
              ["KCl","mM",2500,50,0],
              ["ZnCl2","mM",2.5,0.2,0]]
    DilutionUtil.PrintSolutionSteps(Stats,500,"mL",
                                    BufferName="DI H20")
    print("Pre-treatment dilution")
    Stats = [ ["NiCl2","mM",50,10,0]]
    DilutionUtil.PrintSolutionSteps(Stats,50,"mL",
                                    BufferName="DI H20")
    # roughly speaking, these are (very high, should be able to image in liquid,
    # shouldnt be able to image in liquid), corresponding to...
    # (Pastre, 2003), ns2~0.972, ns2~0.962, ns2~0.953
    MgCl2_concentrations_mM = [120,64,40]
    max_conc = max(MgCl2_concentrations_mM)
    Stats = [ ["2x PRC2 buffer","x",2,1,0],
              ["MgCl2","mM",1000,max(MgCl2_concentrations_mM),0]]
    print("==== MgCl2 Stock ====")
    DilutionUtil.PrintSolutionSteps(Stats,75,"mL",
                                    BufferName="DI H20")
    # serially dilute to what we need..
    DilutionUtil.PrintSerialSteps(Stock=max_conc,Volumes=20,
                                  Desired=MgCl2_concentrations_mM,
                                  ConcString="mM MgCl2",VolString="mL",
                                  BufferString="1x Buffer")
    # try again with NiCl2 if we dont have any luck?
    NiCl2_concentrations_mM = [1,0.05,0.01,5e-3]
    print("==== NiCl2 Stock ====")
    # serially dilute to what we need..
    DilutionUtil.PrintSerialSteps(Stock=10,Volumes=20,
                                  Desired=NiCl2_concentrations_mM,
                                  ConcString="mM NiCl2",VolString="mL",
                                  BufferString="1x Buffer")




if __name__ == "__main__":
    run()
