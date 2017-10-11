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
from GeneralUtil.python import PlotUtilities,CheckpointUtilities,GenUtilities
from GeneralUtil.python.Plot import Scalebar,Annotations,Inset
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

def reload_filtered_landscapes(force_re_filter):
    """
    Returns: tuple of (output of filter_all_landscapes), bin arrays, landscapes
    
    Note: output of filtered landscapes[i] is (bins in filtering i, landscapes 
    in i)
    """
    cache_dir = "./cache/"
    GenUtilities.ensureDirExists(cache_dir)
    n_files = len(GenUtilities.getAllFiles(cache_dir,ext=".pkl"))
    data_dir = "../../EnergyLandscapes/Full_(no_adhesion).pkl_k_31.57pN_nm/"
    # only re-load if we have to 
    if (n_files == 0 or force_re_filter):
        limit = None
    else:
        # just get the key / first landscsape 
        limit = 1
    landscapes = CheckpointUtilities.\
            multi_load(cache_dir=data_dir,load_func=None,force=False,
                       limit=limit)        
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
    return bins,bins_sizes_n,n,list_filtered_n,landscapes
    
class plot_info:
    def __init__(self,key,bin_sizes,average_stdev_energy,stdev_stdev_energy,
                 list_filtered_n):
        self.bin_sizes = bin_sizes
        self.key=key
        self.average_stdev_energy=average_stdev_energy
        self.stdev_stdev_energy=stdev_stdev_energy
        beta = self.key.beta
        kT_to_kcal_per_mol = 0.592
        to_plot_y = lambda x: x *  (beta**2) * 1e9 * (kT_to_kcal_per_mol)**2
        self.average_error_per_bin_plot = to_plot_y(self.average_stdev_energy)
        self.stdev_stdev_energy_per_bin_plot = \
            to_plot_y(self.stdev_stdev_energy)
        self.upper_bound = self.average_error_per_bin_plot+\
                           self.stdev_stdev_energy_per_bin_plot
        self.bin_sizes_nm = bin_sizes * 1e9
        self.min_e = min(self.upper_bound)
        well_max =  1.3*self.min_e
        where_close = np.where(self.upper_bound <= well_max)[0]
        self.idx_n = int(np.round(np.mean(where_close)))
        self.key_filtered = list_filtered_n[self.idx_n][0]
        res_m = bin_sizes[self.idx_n]
        self.res_nm = res_m * 1e9
        self.idx_chosen = np.argmin(np.abs(bin_sizes - res_m))


def get_plotting_info():
    force_re_filter = False
    bins,bin_sizes_n,n,list_filtered_n,landscapes = \
        reload_filtered_landscapes(force_re_filter)
    np.testing.assert_allclose(sorted(bin_sizes_n),n)
    # POST: data is OK, filtered as we want. 
    bin_sizes = np.array([b[1]-b[0] for b in bins])
    max_q = np.max([max(b) for b in bins])
    interp_points = max(n*10)
    # interpolate everything onto a common grid;
    # PRE: all landscapes are offset before filtering.
    x = np.linspace(0,max_q,interp_points)
    # the error in tje
    error_value = lambda f : np.gradient(f.spline_fit.y(x)/(x[1]-x[0]))
    energy_stdev = [ np.std([error_value(f) for f in list_n],axis=0)
                    for list_n in list_filtered_n]
    residuals = [ np.array([l.spline_fit.spline.residual
                            for l in list_n])
                  for list_n in list_filtered_n]
    average_stdev_energy = np.array([np.mean(e)*np.mean(r)
                                     for e,r in zip(energy_stdev,residuals)])
    stdev_stdev_energy = np.array([np.std(e)*np.mean(r)
                                   for e,r in zip(energy_stdev,residuals)])
    key = landscapes[0]
    return plot_info(key=key,average_stdev_energy=average_stdev_energy,
                     stdev_stdev_energy=stdev_stdev_energy,bin_sizes=bin_sizes,
                     list_filtered_n=list_filtered_n)

def run():
    """
    """
    # load all the landscapes here
    inf = CheckpointUtilities.getCheckpoint("./data.pkl",
                                            get_plotting_info,False)
    key = inf.key
    to_x = lambda x: x*1e9
    to_y = lambda y : y * key.beta
    kbT_text = "$k_\mathrm{b}T$"
    eq_str = ("Error =                                             \n"+
              r"$\left(<\sum_{q}(G_q-G_{q,\mathrm{fit}})^2>_N "
              r"\cdot " + \
              r"<\sigma_\frac{dG}{dq}>_N\right)$")
    kbT_text_paren = "(" + kbT_text + ")"
    fig = PlotUtilities.figure(figsize=(3.25,3))
    ax_error = plt.subplot(1,1,1)
    ax_error.set_xscale('log')
    ax_error.set_yscale('log')
    y_plot = inf.average_error_per_bin_plot
    marker_props = dict(color='b',marker='o',mfc="b",markersize=2)
    errorbar_dict = dict(linewidth=0.75,capsize=1.5,
                         **marker_props)
    plt.errorbar(x=inf.bin_sizes_nm,y=y_plot,yerr=0,**errorbar_dict)
    Annotations.relative_annotate(ax=ax_error,s=eq_str,fontsize=8.5,
                                  xy=(0.385,0.8),color='b')
    # plot the called out one 
    chosen_dict = dict(**errorbar_dict)
    chosen_dict['color']='r'
    chosen_dict['mfc']='r'
    chosen_dict['markersize'] *= 1.5
    idx_chosen = inf.idx_chosen
    plt.errorbar(x=inf.bin_sizes_nm[idx_chosen],
                 y=inf.average_error_per_bin_plot[idx_chosen],
                 **chosen_dict)    
    # give one of them a proper error bar
    idx_error = idx_chosen - 16
    plt.errorbar(x=inf.bin_sizes_nm[idx_error],y=y_plot[idx_error],
                 yerr=inf.stdev_stdev_energy_per_bin_plot[idx_error],
                 **errorbar_dict)
    xlim = [1e-2,None]
    error_paren = "((kcal/mol)$^3$/nm)"
    PlotUtilities.lazyLabel("Bin size (nm)",
                            r"Error " + error_paren,"",
                            legend_kwargs=dict(color='b',fontsize=71))
    Annotations.relative_annotate(ax=ax_error,s="{:.1f} nm".format(inf.res_nm),
                                  xy=(inf.res_nm,min(y_plot)*1.2),
                                  color='r',size=6,
                                  xycoords='data')
    plt.xlim(xlim)
    PlotUtilities.savefig(fig,"./filtering.png",pad_inches=0.02)

if __name__ == "__main__":
    run()
