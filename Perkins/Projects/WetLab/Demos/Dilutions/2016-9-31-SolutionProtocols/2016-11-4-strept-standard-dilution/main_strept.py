# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../")
from Util import DilutionUtil 

def run():
    """
    For aliquotting things...
    """
    # TCEP is already present,
    # assume effectively that we want the full aliquot
    # Stats list is formattd like <name,Concentraiton string, stock, desired,
    # already present> s
    ProteinStock = 1
    ProteinDesired = 0.2
    Stats = [ ["TCEP","mM",50,1,0],
              ["Strept (in Aliquot)","mg/mL",ProteinStock,ProteinDesired,0]]
    # get the stocks, desired concntrations, and already-present concentraitons
    Stocks = [s[2] for s in Stats]
    Desired = [s[3] for s in Stats]
    Already = [s[4] for s in Stats]
    Volume = 150 # uL
    Volumes = DilutionUtil.\
        GetVolumesNeededByConcentration(Stocks,Desired,Volume,
                                        AlreadyHaveMass=Already)
    BufferVolume = Volume - sum(Volumes) 
    print("In a total solution of {:.1f}uL...".format(Volume))
    vol_units = "uL"
    for (name,conc_units,conc_stock,desired_conc,_),vol_stock in\
        zip(Stats,Volumes):
        print("\t{:.2f}{:s} of {:.2f}{:s} {:s} for {:.2f}{:s} in solution".\
              format(vol_stock,vol_units,conc_stock,conc_units,name,
                     desired_conc,conc_units))
    print("\tRemainder is ({:.1f}uL) of buffer".format(BufferVolume))
    # Now assume that the actual concentration of Protein was somewhat higher;
    # what do we have in the solution?
    ActualConcentration = 1.9
    FactorMore  = ActualConcentration/ProteinStock
    TrueConcentration = FactorMore * ProteinDesired
    print("For 0.2mg/mL according to nanodrop " + \
          "(assume stock at {:.1f} mg/mL)".format(ActualConcentration))
    DilutionUtil.PrintDilutions([TrueConcentration],[20],[ProteinDesired],
                                UnitConc=["mg/mL"])

                                      

if __name__ == "__main__":
    run()
