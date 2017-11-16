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
    # stock concentration of TCEP, mM
    Stock = 500
    # Desired concentrations
    Desired = [50]
    # desired volumes (for each)
    Volumes = [200]
    buffer_name = "PBS, pH 6.75"
    DilutionUtil.PrintSerialSteps(Stock,Volumes,sorted(Desired)[::-1],
                                  ConcString="mM",BufferString=buffer_name)
    # TCEP is already present,
    # assume effectively that we want the full aliquot
    # Stats list is formattd like <name,Concentraiton string, stock, desired,
    # already present> s
    ProteinStock = 1
    ProteinDesired = 0.2
    Stats = [ ["TCEP","mM",50,1,0],
              ["Strept (in Aliquot)","mg/mL",ProteinStock,ProteinDesired,0]]
    aliquot_volume = 150
    vol_units = "uL"
    n_aliquots_at_a_time = 10
    post_volume = 60
    print("{:d}x aliquots...".format(n_aliquots_at_a_time))
    DilutionUtil.PrintSolutionSteps(Stats,aliquot_volume*n_aliquots_at_a_time,
                                    vol_units,BufferName=buffer_name,
                                    PostVolume=post_volume*n_aliquots_at_a_time)
    print("Single aliquot...")
    DilutionUtil.PrintSolutionSteps(Stats,aliquot_volume,vol_units,
                                    BufferName=buffer_name,
                                    PostVolume=post_volume)
    print("======> Add {:d}uL of {:s} after thawing! <=====".\
          format(post_volume,buffer_name))
                                      

if __name__ == "__main__":
    run()
