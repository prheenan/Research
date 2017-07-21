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
    print("No NiCl2 buffer creation (high salt)")
    Stats = [ ["Tris-HCl","mM",100,10,0],
              ["KCl","mM",1000,300,0],
              ["ZnCl2","mM",2.5,0.1,0]]
    Volume = 250
    vol_units = "mL"
    DilutionUtil.PrintSolutionSteps(Stats,Volume,vol_units,
                                    BufferName="DI H20")
    print("No NiCl2 buffer creation (low salt)")
    Stats = [ ["Tris-HCl","mM",100,10,0],
              ["KCl","mM",1000,25,0],
              ["ZnCl2","mM",2.5,0.1,0]]
    DilutionUtil.PrintSolutionSteps(Stats,Volume,vol_units,
                                    BufferName="DI H20")
    print("Pre-treatment dilution")
    Stats = [ ["Tris-HCl","mM",100,10,0]]
    DilutionUtil.PrintSolutionSteps(Stats,50,vol_units,
                                    BufferName="DI H20")
    # write down how to get the volumes we will want
    volumes = 30
    desired_dilutions_mM = [100,50,25]
    desired_volumes_mL = [volumes for d in desired_dilutions_mM]
    DilutionUtil.PrintSerialSteps(300,desired_volumes_mL,desired_dilutions_mM,
                                  ConcString="mM KCl",VolString="mL",
                                  BufferString="25 mM KCl",
                                  dilution_concentration=25)
    print("=== Low salt dilution === ")
    for desired_dilutions_mM in [0.33,0.083,0.021]:
        desired_volumes_mL = volumes
        DilutionUtil.PrintSerialSteps(3,desired_volumes_mL,
                                      [desired_dilutions_mM],
                                      ConcString="mM NiCl2",VolString="mL",
                                      BufferString="[0mM NiCl2, Same KCl]")



if __name__ == "__main__":
    run()
