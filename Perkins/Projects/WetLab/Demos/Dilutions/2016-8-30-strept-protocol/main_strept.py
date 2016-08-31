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
    stocks = np.array([5,4])
    StreptConcStock = 4.9   # mg/mL
    StreptConcDesired = 0.2 # mg/mL
    VolumeStrept = 10
    VolumeAdd = DilutionUtil.GetVolumeToDilute(StreptConcStock,
                                               VolumeStrept,
                                               StreptConcDesired)
    VolumeTotal = VolumeStrept + VolumeAdd
    TCEPStock = 4   # mM
    TCEPDesired = 1 # mM
    VolumeTCEP = TCEPDesired * VolumeTotal / TCEPStock
    Stats = [ ["TCEP","mM",TCEPStock,VolumeTCEP,"uL",TCEPDesired],
              ["Strept","mg/mL",StreptConcStock,VolumeStrept,"uL",
               StreptConcDesired]]
    BufferVolume = VolumeTotal - VolumeTCEP - VolumeStrept
    print("In a total solution of {:.1f}uL ({:.1f} of buffer)...".\
          format(VolumeTotal,BufferVolume))
    for name,conc_units,conc_stock,vol_stock,vol_units,desired_conc in Stats:
        print("\t{:.1f}{:s} of {:.1f}{:s} {:s} for {:.1f}{:s} in solution".\
              format(vol_stock,vol_units,conc_stock,conc_units,name,
                     desired_conc,conc_units))
    print("\tRemainder {:.1f}uL is buffer".format(BufferVolume))
                                      

if __name__ == "__main__":
    run()
