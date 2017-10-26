# force floating point division. Can still use integer with /
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,copy

sys.path.append("../")
sys.path.append("../../../../../../../../../")
from Research.Perkins.AnalysisUtil.Images import ImageUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

from GeneralUtil.python import CheckpointUtilities,PlotUtilities
from Util import Processing
from skimage.morphology import skeletonize
from skimage.filters import gaussian
from skimage.measure import label,regionprops
from sklearn.metrics.pairwise import pairwise_distances

def is_valid_trace(p):
    coords = p.coords 
    n_coords = coords.shape[0]
    if (n_coords < 2):
        return False
    # POST: something to look at
    expected_shape = (n_coords,n_coords)
    distance_matrix = pairwise_distances(X=coords)
    assert distance_matrix.shape == expected_shape
    lower_is_inf = np.zeros(expected_shape)
    lower_is_inf[np.tril_indices(n=n_coords,k=0)] = np.inf
    upper_idx = np.triu_indices(n=n_coords,k=1)
    lower_is_inf[upper_idx] = distance_matrix[upper_idx]
    # anywherethe distance is >= 2 is also infinite
    max_dist = np.sqrt(2)
    tol = 1e-6
    lower_is_inf[np.where(lower_is_inf >= max_dist + tol)] = np.inf
    best_idx = np.argsort(lower_is_inf,axis=1)
    endpoints = [j for i,row in enumerate(best_idx[:-1])
                 for j in row if lower_is_inf[i,j] < np.inf ]
    edges = [ [i,e] for i,e in enumerate(endpoints)]
    # make sure the last coordinate is represented)
    assert (n_coords -1) in set(endpoints)
    set_of_vertices = set(list([v for list_v in edges for v in list_v]))
    # make sure all the vertices are selected (the first / 0th vector
    # may not be there, since we chose the upper diagonal)
    missing_coords = set(range(1,n_coords))-set(endpoints)
    if len(missing_coords) > 1:
        return False    
    return True

def _filter_single_image(i):
    new_properties = [p
                      for p in i.trace_properties if is_valid_trace(p)]
    return Processing.TraceInfo(raw_image=i.raw_image,
                                trace_properties=new_properties)
    

def filter(images):
    for i_original in images:
        to_ret = _filter_single_image(i_original)
        yield to_ret
        
def run():
    base_dir = Processing.cache_traces()
    cache_dir = Processing.cache_filtered_traces()
    images = CheckpointUtilities.lazy_multi_load(base_dir)
    load_func = lambda: filter(images)
    e = CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_func,
                                       name_func=FEC_Util.fec_name_func,
                                       force=True)
    for tmp in e:
        height_rel = tmp.raw_image.height
        height_rel -= np.median(height_rel)
        fig = PlotUtilities.figure()        
        ax1 = plt.subplot(2,1,1)
        plt.imshow(height_rel,vmax=2e-9,vmin=0)
        PlotUtilities.FormatImageAxis(ax1)
        ax2 = plt.subplot(2,1,2)
        plt.imshow(tmp.label_image,cmap=plt.cm.spectral,interpolation='nearest')
        PlotUtilities.FormatImageAxis(ax2)
        out_name = cache_dir + FEC_Util._fec_name_func(tmp)+".jpeg"
        PlotUtilities.savefig(fig,out_name)
        
if __name__ == "__main__":
    run()

    