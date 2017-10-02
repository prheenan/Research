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
from GeneralUtil.python.Plot import Scalebar,Annotations
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
    for l in landscapes:
        l.offset_to_min()
    res_m = 0.2e-9
    res_nm = res_m * 1e9
    n = sorted(list(set([int(np.ceil(n)) 
                         for n in np.logspace(np.log10(4),4,num=100)])))
    load_func = lambda :  filter_all_landscapes(n,landscapes)
    list_filtered_n = CheckpointUtilities.multi_load("./cache/",load_func,
                                                     force=force_re_filter)
    max_q = max([max(l.q) for l in landscapes])
    bin_sizes = np.array([max_q/n_tmp for n_tmp in n])
    energy_stdev = [ np.std([f.G_0 for f in list_n],axis=0)
                    for list_n in list_filtered_n]
    average_stdev_energy = np.array([np.mean(e) for e in energy_stdev])
    stdev_stdev_energy = np.array([np.std(e) for e in energy_stdev])
    idx_chosen = np.argmin(np.abs(bin_sizes - res_m))
    key = landscapes[0]
    key_filtered = list_filtered_n[idx_chosen][0]
    average_error_per_bin_kT = average_stdev_energy * key.beta / n
    stdev_stdev_energy_per_bin_kT = stdev_stdev_energy * key.beta / n
    bin_sizes_nm = bin_sizes * 1e9
    to_x = lambda x: x*1e9
    to_y = lambda y : y * key.beta
    thresh = 0.05
    kbT_text = "$k_\mathrm{b}T$"
    kbT_text_paren = "(" + kbT_text + ")"
    fig = PlotUtilities.figure()
    ax_error = plt.subplot(2,1,1)
    ax_error.set_xscale('log')
    ax_error.set_yscale('log')
    plt.axhline(thresh,linewidth=1,linestyle='--',color='k',
                label=("{:.2f} ".format(thresh) + kbT_text))
    marker_props = dict(markeredgewidth=0.2,color='b',marker='o',mfc="w",
                        markersize=1.5)
    errorbar_dict = dict(linewidth=0.3,capsize=0.75,elinewidth=0.4,
                         **marker_props)
    plt.errorbar(x=bin_sizes_nm,y=average_error_per_bin_kT,
                 yerr=stdev_stdev_energy_per_bin_kT,**errorbar_dict)
    # plot the called out one 
    chosen_dict = dict(**errorbar_dict)
    chosen_dict['color']='r'
    plt.errorbar(x=bin_sizes_nm[idx_chosen],
                 y=average_error_per_bin_kT[idx_chosen],
                 yerr=stdev_stdev_energy_per_bin_kT[idx_chosen],
                 **chosen_dict)    
    PlotUtilities.lazyLabel("Bin size (nm)",
                            r"<Error per bin> " + kbT_text_paren,"")
    Annotations.relative_annotate(ax=ax_error,s="0.2nm",xy=(0.1,0.1),
                                  color='r',
                                  xycoords='data')
    plt.subplot(2,1,2)
    plt.plot(to_x(key.q),to_y(key.G_0),color='k',alpha=0.3,linewidth=0.5)
    plt.plot(to_x(key_filtered.q),to_y(key_filtered.G_0),color='k',
             linewidth=0.7)
    PlotUtilities.lazyLabel("Extension (nm)",
                            r"$G_0$ " + kbT_text_paren,
                            "Example landscape, filtered to {:.1f} nm".\
                            format(res_nm),
                            title_kwargs=dict(color='r'))
    PlotUtilities.savefig(fig,"./filtering.png")

if __name__ == "__main__":
    run()
