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

def _cache_base(base_dots="../../"):
    return "{:s}Data/".format(base_dots)

def _data_dir(*args,**kw):
    return _cache_base(*args,**kw) + "/Raw/"
    
def cache_raw(*args,**kw):
    return _cache_base(*args,**kw) + "0_cache_raw/"
  
def cache_filtered(*args,**kw):
    return _cache_base(*args,**kw) + "1_cache_filtered/"
   
def cache_retract(*args,**kw): 
    return _cache_base(**kw) + "2_cache_retract/"