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
    _command_line_config

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    data = np.loadtxt('example.csv',delimiter=',',skiprows=3)
    time,sep,force = data[:,0],data[:,1],data[:,2]
    event_indices = _command_line_config.run_feather(in_file='example.csv',
                                                     threshold=1e-2,
                                                     tau=1e-2,
                                                     spring_constant=6.7e-3,
                                                     trigger_time = 0.382,
                                                     dwell_time = 0.992)
    plt.plot(time,force)
    for i in event_indices:
        plt.axvline(time[i])
    plt.show()
    

if __name__ == "__main__":
    run()
