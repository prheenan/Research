# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os,multiprocessing
from shutil import copyfile

sys.path.append("../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util import Learning,Learners,Offline,\
    Scoring,Analysis
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
    cache_directory = "./cache/"
    # limit (per category)
    limit = 200
    n_folds = 5
    pool_size =  multiprocessing.cpu_count()-1
    debugging = False
    copy_files = False
    force_read = False
    force_relearn = False
    force_learn = False
    only_lowest = True
    n_tuning_points = 15
    debug_directory = "./debug_no_event/"
    GenUtilities.ensureDirExists(debug_directory)
    learners_kwargs = dict(n_points_no_event=n_tuning_points,
                           n_points_fovea=n_tuning_points,
                           n_points_wavelet=n_tuning_points)
    positives_directory = InputOutput.get_positives_directory()
    positive_categories = InputOutput.get_categories(positives_directory,
                                                     only_lowest=only_lowest)
    # for each category, predict where events are
    file_name_cache = "{:s}Scores.pkl".format(cache_directory)
    # XXX use just the first N learners
    n_learners = 3
    learners = Learners.get_learners(**learners_kwargs)[:n_learners]
    learners = CheckpointUtilities.\
               getCheckpoint(file_name_cache,Learning.get_cached_folds,
                             force_relearn,positive_categories,
                             force_read,force_learn,cache_directory,limit,
                             n_folds,pool_size=pool_size,
                             learners=learners)
    for l in learners:
        if debugging:
            break
        # XXX determine where things went wrong (load/look at specific examples)
        # plot everything
        best_metric = Offline.best_metric_from_learner(l)
        distance_histogram= Plotting.event_error_kwargs(best_metric)
        Plotting.plot_individual_learner(debug_directory,l,
                                         rupture_kwargs=distance_histogram)
    num_to_plot = limit
    # XXX looking at the worst of the best for the first learner (no event)
    learner = learners[0]
    valid_scores = learner._scores_by_params(train=False)
    best_metric = Offline.best_metric_from_learner(learner)
    best_param_idx = best_metric.best_param_idx
    # get all the scores in the distance ('best case')
    best_x_value = best_metric.x_values[best_param_idx]
    # get the lowest mediancorresponding validation folds
    folds = [f for f in learner.validation_folds[best_param_idx]]
    # get the  folds for the best parameters
    scores = [score for f in folds for score in f.scores]
    # get all the distances
    true_pred = [s.n_true_and_predicted_events() for s in scores]
    median_dist = [s.minimum_distance_median() for s in scores]
    rupture_dist_hists = [s.euclidean_rupture_spectrum_distance()
                          for s in scores]
    good_idx = [i for i,s in enumerate(rupture_dist_hists) if len(s)>0]
    number_relative = [int(abs(t-p)) for t,p in true_pred]
    # get the worst (largest) distances where we arent none
    # XXX note: None is smaller than everything, seems like, so argsort is OK
    sort_idx = good_idx
    # sort from high to low, first elements are most missed and farest off...
    sort_idx = sorted(sort_idx,reverse=True,
                      key=lambda i:(max(rupture_dist_hists[i])))
    worst_n_idx =  sort_idx[:num_to_plot]
    # csv file names are formatted differently 
    debugging_str = ""
    file_names = [scores[i].source_file + debugging_str + scores[i].name 
                  for i in worst_n_idx]
    print([ (number_relative[i],median_dist[i]) for i in worst_n_idx])
    # os.path.split gives <before file,after file>
    load_files = [os.path.basename(f) +".csv.pkl" for f in file_names]
    load_paths_tmp = [(cache_directory + f)
                      for f in load_files]
    # replace the final underscore...
    print("loading: {:s}".format(load_paths_tmp))
    load_paths = []
    for p in load_paths_tmp:
        if (not os.path.isfile(p)):
            print("Couldn't find [{:s}]".format(p))
        else:
            load_paths.append(p)
    examples = [CheckpointUtilities.getCheckpoint(f,None,False) 
                for f in load_paths]
    threshold = best_x_value
    example_numbers = []
    examples_f = [examples[i] for i in example_numbers]
    scores = []
    for i,example in enumerate(examples):
        load_file_name = (os.path.basename(example.Meta.SourceFile) + \
                          example.Meta.Name + ".csv.pkl")
        # copy the pkl file to the debugging location
        load_path = load_file_name
        debugging_file_path = debug_directory + load_file_name
        if (copy_files):
            copyfile(cache_directory + load_file_name,debugging_file_path)
        # get the prediction, save out the plotting information
        example_split,pred_info = \
            Detector._predict_full(example,threshold=threshold)
        score_tmp = Scoring.get_scoring_info(example_split,pred_info.event_idx)
        scores.append(score_tmp)
        # XXX remove...
        continue
        meta = example.Meta
        GenUtilities.ensureDirExists(cache_directory)
        id_data = "{:d}{:s}{:.1f}p={:s}".format(i,meta.Name,meta.Velocity,
                                                str(threshold))
        wave_name = example_split.retract.Meta.Name
        id_string = debug_directory + "db_" + id_data + "_" + wave_name 
        Plotting.debugging_plots(id_string,example_split,pred_info)
    # XXX Debugging
    # get the bcc
    ruptures_valid_true,ruptures_valid_pred = \
        Learning.get_true_and_predicted_ruptures_per_param(learner)  
    true = ruptures_valid_true[best_param_idx]
    pred = ruptures_valid_pred[best_param_idx]
    ruptures_true,loading_true = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(true)
    ruptures_pred,loading_pred = \
        Learning.get_rupture_in_pN_and_loading_in_pN_per_s(pred)    
    _,bins_rupture,_,bins_load = \
        Learning.limits_and_bins_force_and_load(ruptures_pred,ruptures_true,
                                                loading_true,loading_pred,
                                                limit=False)   
    # XXX fix?
    bins_rupture = np.linspace(-10,400,num=15) 
    bins_load = np.linspace(0,np.log10(2e3),num=15)   
    coeffs = Analysis.bc_coeffs_load_force_2d(loading_true,loading_pred,
                                              bins_load,ruptures_true,
                                              ruptures_pred,bins_rupture)
    # just get the 2d (last one
    bcc = 1-coeffs[-1]                                                
    # get the rupture force spectrum stuff
    rupture_dist_hists = [s.euclidean_rupture_spectrum_distance()
                          for s in scores]
    cat_rupture_dist = np.concatenate(rupture_dist_hists)
    cat_dist_safe = cat_rupture_dist[np.where(cat_rupture_dist > 0)]
    PlotUtilities.figure()
    bins = np.logspace(np.log10(min(cat_dist_safe)),
                       np.log10(max(cat_dist_safe)),num=10)
    line = (cat_rupture_dist,bcc,np.percentile(cat_rupture_dist,
                                                  [25,50,75,90,95,97,100]))
    with open("tmp.csv","w") as f:                       
        f.write("\n".join(["{:s}".format(str(l)) for l in line]))
    fig = PlotUtilities.figure()
    plt.hist(cat_rupture_dist,log=True,bins=bins)
    plt.xscale('log')
    PlotUtilities.lazyLabel("RSQ Spectrum Error","Count",
                            "bcc={:.3g}".format(bcc))
    PlotUtilities.savefig(fig,"./out.png")

    # load the worst n back into memory
    # redo the prediction for the worst N, saving to the debug directory



if __name__ == "__main__":
    run()
