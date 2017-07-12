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


sys.path.append("../../../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import CheckpointUtilities,PlotUtilities,GenUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao

import copy 
from matplotlib import gridspec

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = None
    flickering_dir = "../Data/"
    # XXX use the flickering dir for stuff
    cache_dir = flickering_dir 
    GenUtilities.ensureDirExists(flickering_dir)
    force_read_data = False    
    raw_data = IoUtilHao.read_and_cache_data_hao(in_dir,force=force_read_data,
                                                 cache_directory=flickering_dir,
                                                 limit=3)
    example = raw_data[0]                                                 
    example_plot = copy.deepcopy(example)
    example_plot.Force *= 1e12
    example_plot.Separation *= 1e9
    plot_examples = [example_plot]
    x_func = lambda y: y.Separation
    y_func = lambda y: y.Force 
    ylim_pN = [-20,None]
    xlim_nm = [-5,100]
    zoom_regions_nm = [ [23,24],
                        [60,62]]
    # slice the regions 
    regions = [FEC_Util.slice_by_separation(example_plot,*reg) 
               for reg in zoom_regions_nm]
    ylim_pN_zoom = [50,120]
    # # make the plot 
    fig = PlotUtilities.figure((7,3))
    # create the 'top' gridspec
    top_spec = gridspec.GridSpec(1,3)
    # create separate axes for the image and FECs 
    image_spec = gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=top_spec[0])
    data_spec = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=top_spec[1:]) 
    # create a separate axis for the 'zoom in' FEC 
    zoom_in_spec = gridspec.GridSpecFromSubplotSpec(2,len(zoom_regions_nm),
                                                    subplot_spec=data_spec[0]) 
    # finally, create the heat map axis
    heat_map_axis = gridspec.GridSpecFromSubplotSpec(1,1,
                                                     subplot_spec=data_spec[1]) 

    # # plot the image 
    ax = plt.subplot(image_spec[:])
    ax.axis('off')
    xy_text = 0,np.mean(plt.ylim())
    ax.text(*xy_text,s="Cartoon placeholder")
    # # plot the example fec and zoomed regions
    #
    # 'full' example 
    ax_example = plt.subplot(zoom_in_spec[0,:])
    alpha_data = 0.4
    color_data = 'b'
    dict_plot = dict(n_filter_points=2000,
                     style_data=dict(color=color_data,alpha=alpha_data,
                                     linewidth=0.5,linestyle='-'))
    FEC_Plot._fec_base_plot(x_func(example_plot),y_func(example_plot),
                            **dict_plot)
    plt.ylim(ylim_pN)
    plt.xlim(xlim_nm)
    PlotUtilities.lazyLabel("Separation (nm)","Force (pN)","")   
    PlotUtilities.x_label_on_top(ax_example)    
    # pplot all the zoomed regions 
    for i,r in enumerate(regions):
        ax_tmp = plt.subplot(zoom_in_spec[1,i])
        FEC_Plot._fec_base_plot(x_func(r),y_func(r),**dict_plot)
        xlim = plt.xlim()
        PlotUtilities.zoom_effect01(ax_example, ax_tmp,*xlim)
        PlotUtilities.no_x_anything(ax_tmp)
    # # plot the heat map 
    ax3 = plt.subplot(heat_map_axis[:])
    heat_dict = dict(num_bins=(200,800),ConversionOpts=dict())
    FEC_Plot.heat_map_fec(plot_examples,separation_max=max(xlim_nm),
                          **heat_dict)
    plt.ylim(ylim_pN)
    plt.xlim(xlim_nm)                          
    PlotUtilities.lazyLabel("Separation (nm)","Force (pN)","")
    PlotUtilities.savefig(fig,"out.png")

    
if __name__ == "__main__":
    run()
