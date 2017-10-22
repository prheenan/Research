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


class UnfoldingRefolding:
    def __init__(self,unfolding,refolding,info):
        self.unfolding = unfolding
        self.refolding = refolding
        self.info = info
    @property 
    def Meta(self):
        return self.info.Meta

class RefoldingInfo:
    """
    This class encapsulates a refolding experiment, for the purposes of 
    aligning them 
    """
    def __init__(self,data,idx_start,idx_end,spline):
        self.data = data
        self._idx_start_abs= idx_start
        self._idx_end_abs = idx_end
        self.offset = min(idx_start)
        self.idx_start = idx_start - self.offset
        self.idx_end = idx_end - self.offset
        self.spline = spline
    @property        
    def _idx_pairs(self):
        """
        Returns: slice indices marking the start and end of the unfolding
        and refolding events. first is first unfolding pair, second is second
        (refolding) pair, etc
        """
        n_end = len(self.idx_end)
        idx_flat = sorted(list(self.idx_start) + list(self.idx_end))
        slices = [ slice(i,f,1) for i,f in zip(idx_flat[:-1],idx_flat[1:])]
        assert len(slices) == 2*n_end
        assert [s.start for s in slices] + [slices[-1].stop] == idx_flat
        return slices
    def _z_region_property(self,f):
        """
        returns: f, applied to each region of the spline
        """
        slices = self._idx_pairs
        spline_t = self.spline(self.data.Time)
        return [f(spline_t[s]) for s in slices]
    @property 
    def _z_region_max_min(self):
        """
        returns: the maximum over all region minima
        """
        return max(self._z_region_property(min))
    @property 
    def _z_region_min_max(self):
        """
        returns: the minima over all region maxima
        """
        return min(self._z_region_property(max))
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
    return _landscape_base(*args,**kw) + "_0regions/"
    
def cache_landscape_aligned_regions(*args,**kw):
    return _landscape_base(*args,**kw) + "_1regions_aligned/"

    

