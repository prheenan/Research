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

sys.path.append("../../../../../../../../../")
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
from GeneralUtil.python.Plot import Scalebar
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass,WeierstrassUtil
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,mark_inset

def get_filtered_landscapes(bins_array,landscape):
    filtered = []
    for b in bins_array:
        tmp = WeierstrassUtil.\
              _filter_single_landscape(landscape_obj=landscape,bins=b)
        tmp.offset_to_min()
        filtered.append(tmp)
    return filtered

def filter_all_landscapes(n_points_array,landscapes):
    # zero all the landscapes
    max_q = 0
    for l in landscapes:
        l.offset_to_min()
        # get the max
        max_q = max(max_q,max(l.q))
    # get all the filtered landscapes
    bins_array = [np.linspace(0,max_q,num=n) for n in n_points_array]
    filtered_landscapes = [get_filtered_landscapes(bins_array,l)
                           for l in landscapes]
    # transpose: element i is a list of all landscapes of size n_points_array[i]
    to_ret = map(list, zip(*filtered_landscapes))
    return to_ret

def run():
    """
    """
    # load all the landscapes here
    force_re_filter = False
    landscapes = CheckpointUtilities.\
        multi_load(cache_dir="./data/",load_func=None,force=False)
    print(landscapes)
    n = np.logspace(2,4,num=10)
    load_func = lambda :  filter_all_landscapes(n,landscapes)
    list_filtered_n = CheckpointUtilities.multi_load("./cache/",load_func,
                                                     force=force_re_filter)
    energy_stdev = [ np.std([f.G_0 for f in list_filtered_n],axis=1)
                    for list_filtered_n in list_filtered_n]
    fig = PlotUtilities.figure()
    plt.plot(n,energy_stdev)
    PlotUtilities.savefig(fig,"./filtering.png",
                          subplots_adjust=dict(hspace=0.2))

if __name__ == "__main__":
    run()
