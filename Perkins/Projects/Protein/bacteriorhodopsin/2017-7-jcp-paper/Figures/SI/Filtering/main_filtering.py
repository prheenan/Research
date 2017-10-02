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

def get_filtered_landscape(b,landscape):
    tmp = WeierstrassUtil.\
          _filter_single_landscape(landscape_obj=landscape,bins=b)
    tmp.offset_to_min()
    return tmp

def filter_all_landscapes(n_points_array,landscapes):
    # zero all the landscapes
    max_q = 0
    for l in landscapes:
        l.offset_to_min()
        # get the max
        max_q = max(max_q,max(l.q))
    # get all the filtered landscapes
    bins_array = [np.linspace(0,max_q,num=n) for n in n_points_array]
    for b in bins_array:
        yield b,[get_filtered_landscape(b,l) for l in landscapes]            

def run():
    """
    """
    # load all the landscapes here
    force_re_filter = False
    landscapes = CheckpointUtilities.\
        multi_load(cache_dir="../../Energylandscapes/Full_(no_adhesion).pkl/",
                   load_func=None,force=False,limit=10)
    for l in landscapes:
        l.offset_to_min()
    n = sorted(list(set([int(np.ceil(n)) 
                         for n in np.logspace(np.log10(2),4,num=50)])))
    load_func = lambda :  filter_all_landscapes(n,landscapes)
    ret = CheckpointUtilities.multi_load("./cache/",load_func,
                                         force=force_re_filter)
    bins = [r[0] for r in ret]
    list_filtered_n = [r[1] for r in ret]
    bin_sizes_n = [b.size for b in bins]
    sort_idx = np.argsort(bin_sizes_n)
    bins_sizes_n = [bin_sizes_n[i] for i in sort_idx]
    bins = [bins[i] for i in sort_idx]
    list_filtered_n = [list_filtered_n[i] for i in sort_idx]
    np.testing.assert_allclose(sorted(bin_sizes_n),n)
    # POST: data is OK, filtered as we want. 
    max_q = max([max(l.q) for l in landscapes])
    bin_sizes = np.array([max_q/n_tmp for n_tmp in n])
    error_value = lambda f : np.gradient(f.G_0)/np.gradient(f.q)
    energy_stdev = [ np.std([error_value(f) for f in list_n],axis=0)
                    for list_n in list_filtered_n]
    residuals = [ np.array([l.spline_fit_residual
                            for l in list_n])
                  for list_n in list_filtered_n]
    average_stdev_energy = np.array([np.mean(e)*np.mean(r)
                                     for e,r in zip(energy_stdev,residuals)])
    stdev_stdev_energy = np.array([np.std(e)*np.mean(r)
                                   for e,r in zip(energy_stdev,residuals)])
    key = landscapes[0]
    beta = key.beta
    to_plot_y = lambda x: x *  (beta**2) * 1e9
    average_error_per_bin_plot = to_plot_y(average_stdev_energy)
    stdev_stdev_energy_per_bin_plot = to_plot_y(stdev_stdev_energy)
    error_plot_95_pct = stdev_stdev_energy_per_bin_plot*2
    upper_bound = average_error_per_bin_plot+error_plot_95_pct
    bin_sizes_nm = bin_sizes * 1e9
    min_e = min(upper_bound)
    well_max =  1.1*min_e
    where_close = np.where(upper_bound <= well_max)[0]
    idx_n = int(np.round(np.mean(where_close)))
    key_filtered = list_filtered_n[idx_n][0]
    res_m = bin_sizes[idx_n]
    res_nm = res_m * 1e9
    idx_chosen = np.argmin(np.abs(bin_sizes - res_m))
    to_x = lambda x: x*1e9
    to_y = lambda y : y * key.beta
    kbT_text = "$k_\mathrm{b}T$"
    kbT_text_paren = "(" + kbT_text + ")"
    fig = PlotUtilities.figure()
    ax_error = plt.subplot(2,1,1)
    ax_error.set_xscale('log')
    ax_error.set_yscale('log')
    marker_props = dict(markeredgewidth=0.2,color='b',marker='o',mfc="w",
                        markersize=1.5)
    errorbar_dict = dict(linewidth=0.3,capsize=0.75,elinewidth=0.4,
                         lolims=True,**marker_props)
    plt.errorbar(x=bin_sizes_nm,y=average_error_per_bin_plot,
                 yerr=error_plot_95_pct,**errorbar_dict)
    # plot the called out one 
    chosen_dict = dict(**errorbar_dict)
    chosen_dict['color']='r'
    chosen_dict['mfc']='r'
    chosen_dict['markersize'] *= 2
    plt.errorbar(x=bin_sizes_nm[idx_chosen],
                 y=average_error_per_bin_plot[idx_chosen],
                 yerr=error_plot_95_pct[idx_chosen],
                 **chosen_dict)    
    error_paren = "(" + kbT_text+"$^2$/nm)"
    PlotUtilities.lazyLabel("Bin size (nm)",
                            r"Bin error " + error_paren,"")
    Annotations.relative_annotate(ax=ax_error,s="{:.1f}nm".format(res_nm),
                                  xy=(res_nm,well_max*1.2),
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
