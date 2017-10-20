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


class RefoldingInfo:
    """
    This class encapsulates a refolding experiment, for the purposes of 
    aligning them 
    """
    def __init__(self,data,idx_start,idx_end,spline):
        self.data = data
        self.idx_start = idx_start
        self.idx_end = idx_end
        self.spline = spline
    @property
    def Meta(self):
        return self.data.Meta
    

def _cache_base(base_dots="../../"):
    return "{:s}Data/".format(base_dots)
    
def _landscape_base(*args,**kwargs):
    return _cache_base(*args,**kwargs) + "landscape_"
    
def cache_aligned(*args,**kw):
    return _cache_base(*args,**kw) + "4_cache_aligned/"
    
def cache_landscape_regions(*args,**kw):
    return _landscape_base(*args,**kw) + "regions/"

