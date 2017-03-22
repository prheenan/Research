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
sys.path.append("../../../../../")
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Personal.EventDetection.Util import Plotting,InputOutput
from Research.Personal.EventDetection._2SplineEventDetector import Detector

class test_info:
    def __init__(self,file_name):
        self.file_name = file_name + "Concat.csv.pkl"

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = "./"
    data_base = base + "data/"
    debug_directory = "./out/"
    GenUtilities.ensureDirExists(debug_directory)    
    load_paths = GenUtilities.getAllFiles(data_base,ext=".pkl")
    threshold = 0.2
    for i,f in enumerate(load_paths):
        example = CheckpointUtilities.getCheckpoint(f,None,False) 
        # get the prediction, save out the plotting information
        example_split,pred_info = \
            Detector._predict_full(example,threshold=threshold)
        n_found = len(pred_info.event_slices)
        meta = example.Meta
        id_data = "{:d}{:s}{:.1f}p={:s}".format(i,meta.Name,meta.Velocity,
                                                str(threshold))
        wave_name = example_split.retract.Meta.Name
        id_string = debug_directory + "db_" + id_data + "_" + wave_name 
        Plotting.debugging_plots(id_string,example_split,pred_info)
        # after plotting, check if anything went wrong
        n_expected = len(example.Events)
        err_str = "for {:s}, expected {:d}, got {:d}".\
            format(f,n_expected,n_found)
        assert n_found == n_expected , err_str        

if __name__ == "__main__":
    run()
