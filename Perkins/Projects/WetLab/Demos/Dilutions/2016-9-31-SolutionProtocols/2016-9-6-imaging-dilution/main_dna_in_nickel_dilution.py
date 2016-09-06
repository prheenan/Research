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
    # TCEP is already present (25uL at 4mM),
    # assume effectively that we want the full aliquot
    # Stats list is formattd like <name,Concentraiton string, stock, desired,
    # already present> s
    vol_units = "uL"
    Volume = 100 # units of vol_units
    DesiredDNAConcentration = 2
    StockDNAConcentration = 50
    # note: DNA is assumed in 1mM EDTA, so we need to find what volume this
    # translates into; this corresponds to how much extra EDTA we need
    DNAStockVolNeeded = (DesiredDNAConcentration*Volume)/StockDNAConcentration
    Stats = [ ["Ni2+","mM",3,1,-DNAStockVolNeeded],
              ["DNA","ng/uL",StockDNAConcentration,DesiredDNAConcentration,0]]
    DilutionUtil.PrintSolutionSteps(Stats,Volume,vol_units)
    # now do the serial dilutions
    print("For the serial dilutions from {:.1f}{:s}...".\
          format(DesiredDNAConcentration,vol_units))
    DilutionUtil.PrintSerialSteps(Stock=DesiredDNAConcentration,
                                  Volumes=10,
                                  Desired=[2,1,0.5,0.2])

                                      

if __name__ == "__main__":
    run()
