# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../")
from Util import DilutionUtil 

def run():
    """

    """
    pass
    Stats = [ ["TCEP","mM",4,1],
              ["Strept","mg/mL",4.9,0.2],
              ["Pre-aliquotted","Au",2,1]]
    Desired = [s[-1] for s in Stats]
    Stocks = [s[-2] for s in Stats]
    Volume = 200
    Volumes = DilutionUtil.GetVolumesNeededByConcentration(Stocks,
                                                           Desired,Volume)
    BufferVolume = Volume - sum(Volumes)
    print("In a total solution of {:.1f}uL ({:.1f} of buffer)...".\
          format(Volume,BufferVolume))
    vol_units = "uL"
    for (name,conc_units,conc_stock,desired_conc),vol_stock in\
        zip(Stats,Volumes):
        print("\t{:.1f}{:s} of {:.1f}{:s} {:s} for {:.1f}{:s} in solution".\
              format(vol_stock,vol_units,conc_stock,conc_units,name,
                     desired_conc,conc_units))
    print("\tRemainder ({:.1f}uL) is buffer".format(BufferVolume))
                                      

if __name__ == "__main__":
    run()
