# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../..")
from Util import DilutionUtil 

def run():
    """
    For aliquotting things...
    """
    # TCEP is already present (25uL at 4mM),
    # assume effectively that we want the full aliquot
    # Stats list is formattd like <name,Concentraiton string, stock, desired,
    # already present> s
    TotalDNADesiredNanograms = 10e3
    Volume = 136 # units of vol_units
    Stats = [ ["DNA","ng",88.2,TotalDNADesiredNanograms/Volume,0],
              ["loading buffer","x",6,1,0]]
    # get the stocks, desired concntrations, and already-present concentraitons
    Stocks = [s[2] for s in Stats]
    Desired = [s[3] for s in Stats]
    Already = [s[4] for s in Stats]
    vol_units = "uL"
    Volumes = DilutionUtil.\
        GetVolumesNeededByConcentration(Stocks,Desired,Volume,
                                        AlreadyHaveMass=Already)
    BufferVolume = Volume - sum(Volumes) 
    print("In a total solution of {:.1f}uL...".format(Volume))
    for (name,conc_units,conc_stock,desired_conc,_),vol_stock in\
        zip(Stats,Volumes):
        print("\t{:.3g}{:s} of {:.3g}{:s} {:s} for {:.3g}{:s} in solution".\
              format(vol_stock,vol_units,conc_stock,conc_units,name,
                     desired_conc,conc_units))
    print("\tRemainder is ({:.3g}{:s}) of buffer".\
          format(BufferVolume,vol_units))
                                      

if __name__ == "__main__":
    run()
