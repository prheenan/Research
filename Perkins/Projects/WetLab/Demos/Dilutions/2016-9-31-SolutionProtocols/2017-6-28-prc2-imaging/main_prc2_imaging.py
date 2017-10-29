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
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    # per richard, 2017-6-2
    TCEP_final_mM = 25
    print("=== TCEP dilution === ")
    DilutionUtil.PrintSerialSteps(500,[2000],[TCEP_final_mM],
                                  ConcString="x",
                                  BufferString="DI")
    print(" === FLAG buffer ===")
    Stats = [ ["NaCl","M",2000,150,0],
              ["Tris-HCl","mM",100,10,0],
              ["TCEP","mM",25,1,0]]
    Volume = 50
    vol_units = "mL"
    DilutionUtil.PrintSolutionSteps(Stats,Volume,vol_units,
                                    BufferName="DI")
    print("=== Binding buffer dilution === ")
    DilutionUtil.PrintSerialSteps(4,[200],[1],
                                  ConcString="x",
                                  BufferString="DI")
    print("==== DNA (2.5kbp) dilutions ===")
    ConcString = "nM"
    VolString = "uL"
    # stock concentration
    Stock = 1300
    dna_start_conc = 50
    # Desired concentrations
    Desired = [dna_start_conc]
    # desired volumes (for each)
    Volumes = [50]
    DilutionUtil.PrintSerialSteps(Stock,Volumes,Desired,
                                  ConcString=ConcString,
                                  BufferString="TE") 
    print("==== PRC2 dilutions ===")
    ConcString = "uM"
    VolString = "uL"
    # stock concentration
    Stock = 1
    # Desired concentrations for dna
    Volume_total_inc = 8
    dna_desired_conc_nM = 20
    dna_desired_pmol = dna_desired_conc_nM*Volume_total_inc
    ratio_prc2_dilution = 4
    dna_imaging_stock_conc_nM = \
        dna_desired_pmol/(Volume_total_inc)
    # base the PRC2 on it
    Desired_prc2 = ratio_prc2_dilution * dna_imaging_stock_conc_nM/1000 * \
                   np.array([12,6,3,1])
    num_extra = 1
    # desired volumes (for each)
    volume_1x = Volume_total_inc/ratio_prc2_dilution 
    Volumes = [2*volume_1x] + [volume_1x for _ in Desired_prc2[:-1]]
    DilutionUtil.PrintSerialSteps(Stock,Volumes,Desired_prc2,
                                  ConcString=ConcString,
                                  BufferString="1x PRC2") 
    print("=== Binding dilution (note: PRC2 already has 1x binding)===")
    dna_start_conc_inc = 50
    Stats = [ ["PRC2","x",ratio_prc2_dilution,1,0],
              ["4x binding","x",4,1,1*volume_1x],
              ["DNA","nM",dna_start_conc_inc,dna_imaging_stock_conc_nM,0]]
    vol_units = "uL"
    DilutionUtil.PrintSolutionSteps(Stats,Volume_total_inc,vol_units,
                                    BufferName="1x DI")
    # prc2 is in uM in 'Desired_prc2'.
    prc2_final_conc_nM = 1e3 * np.array(Desired_prc2)/ratio_prc2_dilution
    print("(Note: final molar ratios are: {:s})".\
          format(["{:.2f}".format(s) 
                  for s in (1/dna_imaging_stock_conc_nM) * prc2_final_conc_nM]))
    print("Controls: "+
          "\n\t(1)replace Protein by equal volume 1x Buffer")
    print("=== Imaging Dilution ===")
    ConcString = "nM"
    VolString = "uL"
    # stock concentration
    Stock = dna_imaging_stock_conc_nM
    # Desired concentrations
    Desired = [1,0.5,0.2]
    # desired volumes (for each)
    Volumes = [50 for d in Desired]
    DilutionUtil.PrintSerialSteps(Stock,Volumes,Desired,
                                  ConcString=ConcString,
                                  BufferString="NiCl2+ Buf")

if __name__ == "__main__":
    run()
