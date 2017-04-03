# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os
from scipy import interpolate
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from scipy.stats import norm
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities


class ForceExtensionCategory:
    def __init__(self,number,directory=None,sample=None,velocity_nm_s=None,
                 has_events=False,downsample=None):
        self.category_number = number
        self.directory = directory  
        self.velocity_nm_s = velocity_nm_s
        self.sample = sample
        self.has_events = has_events
        self.data = None
        self.downsample_factor =downsample
        self.scores = None
        if (downsample is not None):
            self.is_simulated=True
        else:
            self.is_simulated = False
    def set_scores(self,scores):
        self.scores = scores
    def set_data(self,data):
        """
        sets the pointer to the list of TimeSepForce objects for this category
        
        Args:
            data: list of TimeSepForce objects
        Returns:
            nothing
        """
        self.data = data 



def get_positives_directory(data_base=None):
    """
    reads the (csv) file at file_path, cachine it to cache_directory,
    reading in events 
    
    Args;
        data_base: where the data lives
    Returns
        positive categories base directory
    """
    if (data_base is None):
        data_base = FEC_Util.default_data_root()
    network = data_base
    base_directory = network + "/4Patrick/CuratedData/Masters_CSCI/"
    positives_directory = base_directory + "Positive/650nm-4x-bio/csv/"
    return positives_directory

def read_and_cache_file(file_path,cache_directory,has_events=False,force=True):
    """
    reads the (csv) file at file_path, cachine it to cache_directory,
    reading in events 
    
    Args;
        file_path: where the file lives
        cache_directory: where to cache the file
        has_events:  if this file has events
        force: if true, force a read
    Returns
        TimeSepForce Object
    """
    file_name = os.path.basename(file_path)
    cache_file = cache_directory + file_name+ ".pkl"
    func_to_call = FEC_Util.read_time_sep_force_from_csv
    return CheckpointUtilities.getCheckpoint(cache_file,func_to_call,force,
                                             file_path,has_events=has_events)
                   
def get_category_data(r_obj,force,cache_directory,limit):
    """
    gets the data for a single category data, caching on a per-data basis
    
    Args:
        r_obj: the category object to use
        others:  see set_and_cache_category_data
    Returns:
        list of time,sep,force objects to use
    """
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
        data_file_tmp = read_and_cache_file(f,**kwargs)
        data_in_category.append(data_file_tmp)
        limit_tmp = limit_tmp - 1
        if (limit_tmp <= 0):
            break
    return data_in_category


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
        data = get_category_data(r_obj,force,cache_directory,limit)
        # set the data in this category
        r_obj.set_data(data_in_category)    



def category_read(category,force,cache_directory,limit,debugging=False):
    """
    Reads in all the data associated with a category

    Args:
        category: ForceExtensionCategory object
        force: if true, force re-reading
        cache_directory: if force is not true, where to re-read from
        limit: maximum number to re-read
    Returns:
        list of TimeSepForce objects
    """
    try:
        return get_category_data(category,force,cache_directory,limit)
    except OSError as e:
        if (not debugging):
            raise(e)
        if (category.category_number != 0):
            return []
        print(e)
        # just read in the files that live here XXX just for debugging
        file_names = GenUtilities.getAllFiles(cache_directory,ext="csv.pkl")
        all_files = [CheckpointUtilities.getCheckpoint(f,None,False)
                     for f in file_names]
        return all_files

def simulated_read(downsample_from,category,limit):
    """
    a function which reads in the first [limit] force extension curves from
    downsample_from, slicing the data by category.downsample_factor (assumed 
    > 1)

    Args:
        categories: list of categories to read in
        force_read,cache_directory,limit): see Learning.get_cached_folds
    """
    n_step = int(category.downsample_factor)
    tol = 1e-9
    assert n_step > 1, "simulation should have downfactor > 1"
    assert abs((n_step-int(n_step))) < tol,\
               "simulation should have integer downfactor"
    data = []
    for l in range(limit):
        tmp = downsample_from.data[l]
        slice_v = slice(0,None,n_step)
        data_tmp = FEC_Util.MakeTimeSepForceFromSlice(tmp,slice_v)
        data.append(data_tmp)
    return data
        
def read_categories(categories,force_read,cache_directory,limit):
    """
    a function to read in a most limit force-extension curves, caching as we go

    Args:
        categories: list of categories to read in
        force_read,cache_directory,limit): see Learning.get_cached_folds
    """
    for c in categories:
        # skip simulated categories initially
        if (c.is_simulated):
            continue
        data_tmp = category_read(c,force_read,cache_directory,limit)
        c.set_data(data_tmp)
    # POST: actual data is set up. go ahead and get any simulated data
    # get the lowest loading rate data to downsample
    loading_rates_effective  = [c.velocity_nm_s 
                                if not c.is_simulated else np.inf
                                for c in categories]
    highest_sampled_idx = np.argmin(loading_rates_effective)
    # use the highest sampled
    highest_sampled_category = categories[highest_sampled_idx]
    for c in categories:
        if (not c.is_simulated):
            continue       
        vel_eff = highest_sampled_category.velocity_nm_s*c.downsample_factor
        c.velocity_nm_s = vel_eff
        file_path = "{:s}_sim_{:.1f}".format(cache_directory,c.velocity_nm_s)
        data_tmp = CheckpointUtilities.\
            getCheckpoint(file_path,simulated_read,force_read,  
                          highest_sampled_category,c,limit)
        c.set_data(data_tmp)
    return categories


    
    
def get_categories(positives_directory,use_simulated=False):
    """
    get all the categories associated with the loading rates we will use

    Args:
        positives_directory: base directory where things live
    Returns:
        list of ForceExtensionCategory
    """
    # tuple of <relative directory,sample,velocity> for FEC with events
    max_load = 1000
    positive_meta = \
    [[positives_directory + "1000-nanometers-per-second/","650nm DNA",max_load]]
     #[positives_directory + "500-nanometers-per-second/","650nm DNA",500], 
     #[positives_directory + "100-nanometers-per-second/","650nm DNA",100]]
    # create objects to represent our data categories
    positive_categories = [ForceExtensionCategory(i,*r,has_events=True) 
                           for i,r in enumerate(positive_meta)]
    if (use_simulated):
        downsample_factors = sorted([2,3,4,10,20,100,1000])
        kw = lambda i: dict(number=(len( positive_categories) + i))
        simulated_categories = [ForceExtensionCategory(downsample=d,**kw(i))
                                for i,d in enumerate(downsample_factors)]
    else:
        simulated_categories = []
    return simulated_categories[::-1] + positive_categories

    
