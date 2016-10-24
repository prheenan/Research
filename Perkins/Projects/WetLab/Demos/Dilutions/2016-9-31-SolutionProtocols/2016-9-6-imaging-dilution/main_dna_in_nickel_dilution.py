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
    # Stats list is formattd like <name,Concentraiton string, stock, desired,
    # already present> 
    vol_units = "uL"
    DesiredDNAConcentration = 2
    # Pre dilution concntration
    PreDilutionConcentration = 68
    # sotkc concentration is what we *want* the stock to be at, for imaging
    StockDNAConcentration = 30
    stocks = np.array([PreDilutionConcentration])
    # what volume are the stocks, in uL
    volumes = np.array([68])
    # The serial aguments
    SerialArgs = dict(Stock=DesiredDNAConcentration,
                      Volumes=20,
                      Desired=[2,0.5])
    Volume = DilutionUtil.StockVolumeNeededForSerialDilution(**SerialArgs)
    # what the post-dilution concenration is, ng/uL
    DesiredConc = np.array([StockDNAConcentration])
    obj = DilutionUtil.PrintDilutions(stocks,volumes,DesiredConc)
                                      
    # note: DNA is assumed in 1mM EDTA, so we need to find what volume this
    # translates into; this corresponds to how much extra EDTA we need. This
    # is effectively a *negative* mass present, (1mM)*ul for each uL 
    DNAStockVolNeeded = (DesiredDNAConcentration*Volume)/StockDNAConcentration
    Stats = [ ["Ni2+","mM",3,1,(-1) * DNAStockVolNeeded],
              ["TE-DNA","ng/uL",StockDNAConcentration,
               DesiredDNAConcentration,0]]
    DilutionUtil.PrintSolutionSteps(Stats,Volume,vol_units,
                                    BufferName="10mM Hepes, pH7")
    # now do the serial dilutions
    print("For serial dilutions from {:.1f}ng/uL into 1mM NiCl2,10mM Hepes...".\
          format(DesiredDNAConcentration))
    DilutionUtil.PrintSerialSteps(**SerialArgs)

                                      

if __name__ == "__main__":
    run()
