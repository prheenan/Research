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
from Research.Personal.EventDetection.Util import Analysis,Plotting

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
            data_in_category.append(read_and_cache_file(f,**kwargs))
            limit_tmp = limit_tmp - 1
            if (limit_tmp == 0):
                break
        # set the data in this category
        r_obj.set_data(data_in_category)    
        
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
    positives_directory = base_directory + "Positive/"
    negatives_directory = "XXX TODO"
    cache_directory = "./cache/"
    # tuple of <relative directory,sample,velocity> for FEC with events
    positive_meta = \
      [[positives_directory + "650nm-4x-bio/csv/500-nanometers-per-second/",
        "650nm DNA",500]]
    # tuple of <relative directory,sample,velocity> for FEC without events
    negative_meta = \
      [[negatives_directory + "/500-nanometers-per-second/csv/",
        "Negative Control",500]]
    # create objects to represent our data categories
    positive_categories = [ForceExtensionCategory(*r,has_events=True) 
                           for r in positive_meta]
    force = False
    # limit (per category)
    limit = 1
    # get the positive events
    set_and_cache_category_data(positive_categories,
                                cache_directory=cache_directory,force=force,
                                limit=limit)
    example = positive_categories[0].data[0]
    example_split = Analysis.split_FEC_by_meta(example)
    approach = example_split.approach
    retract = example_split.retract 
    # get the autocorrelation time of the retract force (what we care about)
    x,f = retract.Time,retract.Force
    dx = np.median(np.diff(x))
    tau,auto_coeffs,auto_correlation = Analysis.auto_correlation_tau(x,f)
    num_points = int(np.round(tau/dx))
    # zero out everything to the approach using the autocorrelation time 
    Analysis.zero_by_approach(example_split,num_points)
    # get an interpolator for the retract force and separation
    force_interpolator = Analysis.spline_interpolator(tau,x,f)
    separation_interpolate = Analysis.spline_interpolator(tau,x,
                                                          retract.Separation)
    # get the residual mean and standard deviation, from the spline...
    f_interp_at_x = force_interpolator(x)
    mu,std = Analysis.spline_residual_mean_and_stdev(f,f_interp_at_x)
    force_cdf = Analysis.spline_gaussian_cdf(f,f_interp_at_x,std)
    force_cdf_complement = 1-force_cdf
    # make a threshold in probability (this will likely be machine-learned) 
    thresh = 1e-4
    # plot everything
    fig = PlotUtilities.figure(figsize=(8,16))
    Plotting.plot_distribution(x,f,f_interp_at_x,force_cdf,thresh,
                               num_bins=500)
    PlotUtilities.savefig(fig,cache_directory + "out.png")
    # get the negative events
    # XXX 


if __name__ == "__main__":
    run()
