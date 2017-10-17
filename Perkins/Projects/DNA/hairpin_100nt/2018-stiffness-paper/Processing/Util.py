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

def _cache_base():
    return "../../Data/"

def _data_dir():
    return _cache_base() + "/Raw/"
    
def cache_raw():
    return _cache_base() + "0_cache_raw/"
  
def cache_filtered():
    return _cache_base() + "1_cache_filtered/"
   
def cache_retract(): 
    return _cache_base() + "2_cache_retract/"