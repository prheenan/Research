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

class TraceInfo:
    def __init__(self,raw_image,trace_properties):
        self.raw_image = raw_image
        self.trace_properties = trace_properties
    @property
    def label_image(self):
        to_ret = np.zeros_like(self.raw_image.height)
        for i,p in enumerate(self.trace_properties):
            for x,y in p.coords:
                to_ret[x,y] = (i+1)
        return to_ret        
    @property
    def Meta(self):
        return self.raw_image.Meta

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
    
def cache_traces(*args,**kw):
    return _cache_dir(2,"TracedImages")
    
def cache_filtered_traces(*args,**kw):
    return _cache_dir(3,"FilteredTraces")        