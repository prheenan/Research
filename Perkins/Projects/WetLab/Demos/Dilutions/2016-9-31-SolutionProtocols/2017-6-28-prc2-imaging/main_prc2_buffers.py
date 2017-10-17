# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,copy


sys.path.append("../../../../")
from Util import DilutionUtil 

def run():
    print("4x buffer creation")
    # see: richard, 2017-6-27. This is the 1x buffer, multiply by 4 below
    """
    Here is the 1X [4x] binding buffer that you should use:
    50    [200] mM Tris-HCl pH 7.5 at 25C
    25    [100] mM KCl
    0.1   [0.4] mM ZnCl2
    2     [8]   mM 2-mercaptoethanol
    0.05% [0.2%]v/v NP-40
    5%    [20%] v/v glycerol
    """
    HEPES_1x = 50
    KCL_1x = 25
    ZnCl2_1x = 0.1
    Stats = [ ["HEPES","mM",975,HEPES_1x,0],
              ["KCl","mM",2500,KCL_1x,0],
              ["ZnCl2","mM",15,ZnCl2_1x,0],
              ["2-mercaptoethanol","mM",201,2,0],
              ["NP-40","%v",10,0.05,0],
              ["Glycerol","%v",100,5,0],
    ]
    # I convert it into the 4x buffer
    conc_mult = 4
    stats_4x = copy.deepcopy(Stats)
    for i in range(len(Stats)):
        stats_4x[i][3] *= conc_mult
    DilutionUtil.PrintSolutionSteps(stats_4x,10,"mL",
                                    BufferName="DI H20")
    print("pH 7.5 buffer, 2x...")
    # make a 2x buffer for AFM imaging with just the salts
    Stats = [ ["HEPES","mM",500,2*HEPES_1x,0],
              ["KCl","mM",2500,2*KCL_1x,0],
              ["ZnCl2","mM",1.25,2*ZnCl2_1x,0]]
    DilutionUtil.PrintSolutionSteps(Stats,500,"mL",
                                    BufferName="DI H20")

    print("No Divalent buffer creation, 2x, HEPES")
    KCl = 50
    Stats = [ ["HEPES","mM",500,20,0],
              ["KCl","mM",2500,KCl,0],
              ["ZnCl2","mM",1.25,0.2,0]]
    DilutionUtil.PrintSolutionSteps(Stats,500,"mL",
                                        BufferName="DI H20")
    print("...")
    print("No Divalent buffer creation, 2x")
    Stats = [ ["Tris-HCl","mM",1000,10,0],
              ["KCl","mM",2500,300,0],
              ["ZnCl2","mM",1.25,0.2,0]]
    DilutionUtil.PrintSolutionSteps(Stats,50,"mL",
                                    BufferName="DI H20")
    print("High salt buffer creation, 1x")
    Stats = [ ["Tris-HCl","mM",1000,10,0],
              ["KCl","mM",2500,500,0],
              ["ZnCl2","mM",1.25,0.1,0]]
    DilutionUtil.PrintSolutionSteps(Stats,50,"mL",
                                    BufferName="DI H20")

    print("...")

    max_conc = 120
    Stats = [ ["2x PRC2 buffer","x",2,1,0],
              ["MgCl2","mM",1000,max_conc,0]]
    DilutionUtil.PrintSolutionSteps(Stats,50,"mL",
                                    BufferName="DI H20")
    # use NiCl2 assay 
    Stats = [ ["2x PRC2 buffer","x",2,1,0],
              ["NiCl2","mM",25,3,0]]
    DilutionUtil.PrintSolutionSteps(Stats,50,"mL",
                                    BufferName="DI H20")
    # serially dilute to what we need..
    # roughly speaking, these are (very high, should be able to image in liquid,
    # shouldnt be able to image in liquid), corresponding to...
    # (Pastre, 2003), ns2~0.972, ns2~0.962, ns2~0.953
    MgCl2_concentrations_mM = [5]
    DilutionUtil.PrintSerialSteps(Stock=max_conc,Volumes=20,
                                  Desired=MgCl2_concentrations_mM,
                                  ConcString="mM MgCl2",VolString="mL",
                                  BufferString="1x Buffer")
    print("==== NiCl2 Stock ====")
    # serially dilute to what we need..
    DilutionUtil.PrintSerialSteps(Stock=25,Volumes=30,
                                  Desired=[0.3],
                                  ConcString="mM NiCl2",VolString="mL",
                                  BufferString="1x Buffer")




if __name__ == "__main__":
    run()
