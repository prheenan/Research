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

def _base_dir(base="../../"):
    return base + "Data/"

def _base_cache_dir(*args,**kw):
    return _base_dir(*args,**kw) + "Cache/"
    
def _cache_dir(num,str,*args,**kw):    
    return _base_cache_dir(*args,**kw) + "_{:d}_{:s}/".format(num,str)
    
def raw_data_dir(*args,**kw):
    return _base_dir(*args,**kw) + "Raw/"
    
def cache_raw_images(*args,**kw):
    return _cache_dir(0,"Images")
    
def cache_select_images(*args,**kw):
    return _cache_dir(1,"SelectedImages")