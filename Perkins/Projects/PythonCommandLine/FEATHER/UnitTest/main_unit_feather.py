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
sys.path.append("../../../../../../")

from Research.Personal.EventDetection._2SplineEventDetector import \
    _command_line_config,Detector

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    # there are a couple of ways to call FEATHER in python. 
    # # (1) through an intermediate '.csv' file
    data = np.loadtxt('example.csv',delimiter=',',skiprows=0)
    time,sep,force = data[:,0],data[:,1],data[:,2]
    threshold = 1e-3
    tau = 2e-2
    meta_dict = dict(threshold=threshold,
                     tau=tau,
                     spring_constant=6.7e-3,
                     trigger_time = 0.382,
                     dwell_time = 0.992)
    event_indices_1 = _command_line_config.run_feather(in_file='example.csv',
                                                       **meta_dict)
    # # (2) directly, using python arrays and a constructed fec object. This is
    # #    likely to be much faster, since there is no extra file IO.
    fec = _command_line_config.make_fec(time=time,separation=sep,force=force,
                                        **meta_dict)
    event_indices_2 = _command_line_config.predict_indices(fec,
                                                           tau_fraction=tau,
                                                           threshold=threshold)
    # make sure the two methods are consistent
    assert np.allclose(event_indices_1,event_indices_2) , \
        "Programming error; FEATHEr methods are identical"
    # POST: they are consistent. go ahead and plot force vs time, add lines
    # where an event is predicted
    plt.plot(time,force*(-1))
    for i in event_indices_1:
        plt.axvline(time[i])
    plt.show()
    print("Found events at indices: {:s}".format(event_indices_1))
    

if __name__ == "__main__":
    run()
