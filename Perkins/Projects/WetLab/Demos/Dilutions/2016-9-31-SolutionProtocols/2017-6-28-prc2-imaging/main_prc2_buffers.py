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
    """
    # see: richard, 2017-6-27. This is the 1x buffer, multiply by 4 below

    Here is the 1X [4x] binding buffer that you should use:
    50    [200] mM Tris-HCl pH 7.5 at 25C
    25    [100] mM KCl
    0.1   [0.4] mM ZnCl2
    2     [8]   mM 2-mercaptoethanol
    0.05% [0.2%]v/v NP-40
    5%    [20%] v/v glycerol

    # from Xueyin:
    #  I think you could eliminate NP40, which initially was put to prevent aggregation. PRC2 is fine without glycerol but I am not sure about DNA-PRC2 interaction. I would imagine it is fine as Glycerol is for keeping the protein stable in most cases.
    # so the new buffer becomes (I substitute TCEP for beta, since the protein
    # is stored in 1mM TCEP). HEPES is better for imaging, and I'm just
    # going to do 50mM to keep everything consistent.

    Here is the 1X [4x] binding buffer that you should use:
    50    [200] mM HEPES pH 7.5 at 25C
    25    [100] mM KCl
    0.1   [0.4] mM ZnCl2
    1     [4]   mM TCEP
    """
    HEPES_1x = 50
    KCL_1x = 25
    ZnCl2_1x = 0.1
    Stats = [ ["HEPES","mM",975,HEPES_1x,0],
              ["KCl","mM",2500,KCL_1x,0],
              ["ZnCl2","mM",15,ZnCl2_1x,0],
              ["TCEP","mM",50,1,0],
    ]
    # I convert it into the 4x buffer
    conc_mult = 4
    stats_4x = copy.deepcopy(Stats)
    for i in range(len(Stats)):
        stats_4x[i][3] *= conc_mult
    DilutionUtil.PrintSolutionSteps(stats_4x,50,"mL",
                                    BufferName="DI H20")
    print("No Divalent buffer creation, 2x, HEPES")
    KCl = 50
    volume_mL = 50
    Stats = [ ["HEPES","mM",500,20,0],
              ["KCl","mM",2500,KCl,0]]
    DilutionUtil.PrintSolutionSteps(Stats,volume_mL,"mL",
                                        BufferName="DI H20")


if __name__ == "__main__":
    run()
