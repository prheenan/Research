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
    # volume for a single aliquot
    VolumeAliquot_uL = 150
    ProteinStock_mg_mL = 0.5
    ProteinDesired_mg_mL = 0.2
    # 1mL total given
    StockProteinVolume_uL = 1e3
    # write down how to dilute the TCEP
    Buffer= "PBS, pH 6.75"
    DesiredTCEPStock_mM = 50
    DesiredTCEPAliquot_mM = 1
    # determine how many aliquots we can make
    NumAliquots = StockProteinVolume_uL/VolumeAliquot_uL
    TotalVolumeTCEP = 50
    DilutionUtil.PrintSerialSteps(500,Volumes=[TotalVolumeTCEP],
                                  Desired=[DesiredTCEPStock_mM],ConcString="mM",
                                  BufferString=Buffer)
    # Stats list is formattd like <name,Concentraiton string, stock, desired,
    # already present>
    Stats = [ ["TCEP","mM",DesiredTCEPStock_mM,DesiredTCEPAliquot_mM,0],
              ["T-Strept","mg/mL",ProteinStock_mg_mL,ProteinDesired_mg_mL,0]]
    DilutionUtil.PrintSolutionSteps(Stats,VolumeAliquot_uL,vol_units="uL",
                                    BufferName=Buffer,PostVolume=60)

                                      

if __name__ == "__main__":
    run()
