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

