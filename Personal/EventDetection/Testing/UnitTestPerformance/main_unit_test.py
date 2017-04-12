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
from Research.Personal.EventDetection.Util import Plotting,InputOutput,Scoring,\
    Learning,Analysis
from Research.Personal.EventDetection._2SplineEventDetector import Detector

def check_bcc(examples,predicted,bcc_threshold=0.0434,
              rupture_tuple=(0.0699,0.446)):
    # get the scoring objects
    scores = []
    for example_split,pred_info in zip(examples,predicted):          
        score = Scoring.get_scoring_info(example_split,pred_info.event_idx)
        scores.append(score)
    # get the true and predicted rupture objects 
    true = [obj for s in scores for obj in s.ruptures_true]
    pred = [obj for s in scores for obj in s.ruptures_predicted]
    ruptures_true,loading_true = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(true)
    ruptures_pred,loading_pred = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(pred)    
    _,bins_rupture,_,bins_load = \
        Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                loading_true,loading_pred,
                                                limit=True)   
    coeffs = Analysis.bc_coeffs_load_force_2d(loading_true,loading_pred,
                                              bins_load,ruptures_true,
                                              ruptures_pred,bins_rupture)
    # just get the 2d (last one
    bcc = 1-coeffs[-1]          
    bcc_str =  "bcc is {:.4g}".format(bcc)
    assert bcc <= bcc_threshold , bcc_str
    print(bcc_str)
    # get the rupture force spectrum stuff
    rupture_dist_hists = [s.euclidean_rupture_spectrum_distance()
                          for s in scores]
    cat_rupture_dist = np.concatenate(rupture_dist_hists)
    median_rupture_dist = np.median(cat_rupture_dist)
    stdev_rupture_dist = np.std(cat_rupture_dist)
    tuple_true = (median_rupture_dist,stdev_rupture_dist)
    tuple_str = "rupture distance/stdev is {:s}".format(tuple_true)
    assert (np.array(tuple_true) <= np.array(rupture_tuple)).all() , tuple_str
    print(tuple_str)

def check_single_file(example_split,pred_info,fractional_error_tolerance):
    n_found = len(pred_info.event_slices)
    meta = example_split.retract.Meta
    # after plotting, check if anything went wrong
    events=  example_split.get_retract_event_centers()
    n_expected = len(events)
    err_str = "for {:s}, expected {:d}, got {:d}".\
        format(meta.Name,n_expected,n_found)
    assert n_found == n_expected , err_str   
    # POST: number of events match. check that the locations match (within
    # fractional_error_tolerance * number of points)
    predicted_centers = np.array(sorted(pred_info.event_idx))
    actual_centers = np.array(sorted(events))
    errors = abs(predicted_centers-actual_centers)
    n = example_split.retract.Force.size
    rel_errors = [_/n for _ in errors]
    err_str = ("Fractional Error(s): " + \
               ",".join(["{:.3g}".format(_) for _ in rel_errors]))
    for e in errors:
        assert e <= fractional_error_tolerance * n , err_str
    print(err_str)
    return rel_errors
               
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
    threshold = 1e-3
    fractional_error_tolerance = 4.04e-3
    predicted,examples = [],[]
    max_error = 0
    for i,f in enumerate(load_paths[1:]):
        example = CheckpointUtilities.getCheckpoint(f,None,False) 
        # get the prediction, save out the plotting information
        example_split,pred_info = \
            Detector._predict_full(example,threshold=threshold)
        # save out a debugging plot (before predicting/maybe failing)
        meta = example_split.retract.Meta
        id_data = "{:d}{:s}{:.1f}p={:s}".format(i,meta.Name,meta.Velocity,
                                                str(threshold))
        wave_name = meta.Name
        id_string = debug_directory + "db_" + id_data + "_" + wave_name 
        Plotting.debugging_plots(id_string,example_split,pred_info)
        kw = dict(fractional_error_tolerance=fractional_error_tolerance)
        errors = check_single_file(example_split,pred_info,**kw)
        max_error = max(max_error,max(errors))
        # save the info so we can check for the bcc
        predicted.append(pred_info)
        examples.append(example_split)
    print("The maximum relative error was {:.4g}".format(max_error))
    # POST: looks okay, but lets just the bcc.
    check_bcc(examples,predicted)

if __name__ == "__main__":
    run()
