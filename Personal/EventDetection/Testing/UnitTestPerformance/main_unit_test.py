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
    threshold = 1e-2
    fractional_error_tolerance = 0.01    
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
        # POST: number of events match. check that the locations match (within
        # fractional_error_tolerance * number of points)
        predicted_centers = np.array(sorted(pred_info.event_idx))
        actual_centers = \
            np.array(sorted(example_split.get_retract_event_centers()))
        errors = abs(predicted_centers-actual_centers)
        n = example_split.retract.Force.size
        for e in errors:
            assert e <= fractional_error_tolerance * n
        print(("Fractional Error(s): " + \
               ",".join(["{:.3g}".format(_/n) for _ in errors])))

if __name__ == "__main__":
    run()
