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

sys.path.append("../../../../../..")
from GeneralUtil.python import PlotUtilities
from IgorUtil.PythonAdapter import PxpLoader

import allantools,copy # https://github.com/aewallin/allantools/

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    file_name = "defl.ibw"
    defl_volts_wave = PxpLoader.read_ibw_as_wave(file_name)
    defl_volts = defl_volts_wave.DataY
    invols_nm_per_volt = 58.75
    spring_pN_per_nM = 6.84
    force_pN = defl_volts * invols_nm_per_volt * spring_pN_per_nM
    # zero out to the maximum
    force_pN -= np.mean(force_pN)
    data_interval_in_s = 1/500e3
    rate = 1/float(data_interval_in_s) # data rate in Hz
    taus = np.logspace(-5,2,num=50,base=10) #  tau-values in seconds
    # fractional frequency data
    (taus_used, adev, adeverror, adev_n) = \
        allantools.adev(force_pN, data_type='freq', rate=rate, taus=taus)
    fig = PlotUtilities.figure()
    plt.errorbar(x=taus_used,y=adev,yerr=adeverror)
    plt.xscale('log')
    plt.yscale('log')
    plt.axhline(1,color='r',linestyle='--',label="1 pN")
    PlotUtilities.lazyLabel("Averaging time (s)","Force (pN)","")
    PlotUtilities.savefig(fig,"out.png")

if __name__ == "__main__":
    run()
