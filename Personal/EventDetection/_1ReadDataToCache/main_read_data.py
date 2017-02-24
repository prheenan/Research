# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os

sys.path.append("../../../../")
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from GeneralUtil.python.IgorUtil import SavitskyFilter
from Research.Personal.EventDetection.Util import Analysis,Plotting,InputOutput
from Research.Personal.EventDetection._2SplineEventDetector import Detector

class ForceExtensionCategory:
    def __init__(self,directory,velocity_nm_s,sample,has_events):
        self.directory = directory  
        self.velocity_nm_s = velocity_nm_s
        self.sample = sample
        self.has_events = has_events
        self.data = None
    def set_data(self,data):
        """
        sets the pointer to the list of TimeSepForce objects for this category
        
        Args:
            data: list of TimeSepForce objects
        Returns:
            nothing
        """
        self.data = data
                          
def set_and_cache_category_data(categories,force,cache_directory,limit):
    """
    loops through each category, reading in at most limit files per category,
    caching the csv files to cache_directory
    
    Args:
        categories: list of ForceExtensionCategory objects. will have their
        data set with the appropriate TimeSepForce objects 
        
        force: whether to force re-reading
        cache_directory: where to save the files
        limit: maximum number of files to each (per category)
    Returns:
        nothing, but sets the data of set_categories
    """
    for i,r_obj in enumerate(categories):
        # restart the limit each category
        limit_tmp = limit
        data_in_category = []
        # get the label for this dataset.
        dir_v = r_obj.directory
        all_files = GenUtilities.getAllFiles(dir_v,ext=".csv")
        kwargs =dict(cache_directory=cache_directory,
                     has_events = r_obj.has_events,force=force)
        # reach all the files until we reach the limit
        for f in all_files:
            data_file_tmp = InputOutput.read_and_cache_file(f,**kwargs)
            data_in_category.append(data_file_tmp)
            limit_tmp = limit_tmp - 1
            if (limit_tmp == 0):
                break
        # set the data in this category
        r_obj.set_data(data_in_category)    
        
def debug_plotting(example,cache_directory,out_file_name):    
    example_split = Analysis.split_FEC_by_meta(example)
    approach = example_split.approach
    retract = example_split.retract 
    # get the autocorrelation time of the retract force (what we care about)
    x,f = retract.Time,retract.Force
    separation = retract.Separation
    dx = np.median(np.diff(x))
    deg_auto = 1
    tau,auto_coeffs,auto_correlation = Analysis.\
        auto_correlation_tau(x,f,deg_autocorrelation=deg_auto)
    num_points = int(np.ceil(tau/dx))
    # zero out everything to the approach using the autocorrelation time 
    Analysis.zero_by_approach(example_split,num_points)
    # XXX only look at after the nominal zero point?
    # get an interpolator for the retract force and separation
    force_interpolator = Analysis.spline_interpolator(tau,x,f)
    separation_interpolate = Analysis.spline_interpolator(tau,x,separation)
    # get the residual mean and standard deviation, from the spline...
    f_interp_at_x = force_interpolator(x)
    mu,std = Analysis.spline_residual_mean_and_stdev(f,f_interp_at_x)
    force_cdf = Analysis.spline_gaussian_cdf(f,f_interp_at_x,std)
    force_cdf_complement = 1-force_cdf
    # get the derivative of the splined data
    derivative_at_x = force_interpolator.derivative()(x)
    # make a threshold in probability (this will likely be machine-learned) 
    thresh = 1e-4
    # plot everything
    fig = PlotUtilities.figure(figsize=(8,20))
    Plotting.plot_autocorrelation_log(x, tau,auto_coeffs,auto_correlation)
    PlotUtilities.savefig(fig,cache_directory + out_file_name + "out.png")    
    # XXX this is a *bad* way of doing things (should pass example_split in)
    return example_split
        
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    network = FEC_Util.default_data_root()
    base_directory = network + "/4Patrick/CuratedData/Masters_CSCI/"
    positives_directory = base_directory + "Positive/650nm-4x-bio/csv/"
    negatives_directory = "XXX TODO"
    cache_directory = "./cache/"
    # tuple of <relative directory,sample,velocity> for FEC with events
    positive_meta = \
      [[positives_directory + "1000-nanometers-per-second/","650nm DNA",100],
       [positives_directory + "500-nanometers-per-second/","650nm DNA",500], 
       [positives_directory + "100-nanometers-per-second/","650nm DNA",1000]]
    # tuple of <relative directory,sample,velocity> for FEC without events
    negative_meta = \
      [[negatives_directory + "/500-nanometers-per-second/csv/",
        "Negative Control",500]]
    # create objects to represent our data categories
    positive_categories = [ForceExtensionCategory(*r,has_events=True) 
                           for r in positive_meta]
    force = False
    # limit (per category)
    limit = 5
    # get the positive events
    set_and_cache_category_data(positive_categories,
                                cache_directory=cache_directory,force=force,
                                limit=limit)
    thresh = 1e-2                                
    # for each category, predict where events are
    for i,category in enumerate(positive_categories):
        for j,example in enumerate(category.data):
            id = "{:d}_{:d}".format(i,j)
            example_split = debug_plotting(example,cache_directory,id)
            m_func =Detector.adhesion_function_for_split_fec(example_split)
            info = Detector._predict_helper(example_split,threshold=thresh,
                                            condition_function=m_func)
            # XXX fix threshhold
            fig = PlotUtilities.figure(figsize=(8,12))    
            Plotting.plot_prediction_info(example_split,info)
            PlotUtilities.savefig(fig,cache_directory + "info" + id + ".png")

    # get the negative events
    # XXX 


if __name__ == "__main__":
    run()
