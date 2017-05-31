# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,h5py 

sys.path.append("../../../../../")

from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util 
from IgorUtil.PythonAdapter import TimeSepForceObj

def read_matlab_file_into_fec(input_file):
    f = h5py.File(input_file,'r') 
    get = lambda x: f[x].value
    # get the FEC data
    time = get('time')
    separation = get('separation')
    force = get('force')
    # get the meta data needed
    dwell_time = get('dwell_time')
    trigger_time = get('trigger_time')
    meta_dict = dict(DwellTime=dwell_time,
                     TriggerTime=trigger_time,
                     DwellSetting=1,
                     # XXX
                     Name=input_file,
                     K=1e-3,
                     # invols is not used, se we just set it to one
                     Invols=1)
    # XXX approach/retract velocity, spring constant 
    data = TimeSepForceObj.data_obj_by_columns_and_dict(time=time,
                                                        sep=separation,
                                                        force=force,
                                                        meta_dict=meta_dict)
    to_ret = TimeSepForceObj.TimeSepForceObj()
    to_ret.LowResData = data
    plt.plot(to_ret.Time[-10000:],to_ret.Force[-10000:],'r,')


def run():
    """
    """
    read_matlab_file_into_fec('tmp.mat')

if __name__ == "__main__":
    run()
