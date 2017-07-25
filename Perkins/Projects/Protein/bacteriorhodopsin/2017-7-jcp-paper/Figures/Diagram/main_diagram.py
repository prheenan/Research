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
from GeneralUtil.python.Plot import Scalebar
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
    zoom_regions_nm = [ [21,23],
                        [59,62]]
    adhesion_max_nm = 19
    # plot the helical regions...
    regions_nm = [ [[adhesion_max_nm,30],"ED Helix",'royalblue'],
                   [[31,45],"CB Helix",'orangered'],
                   [[55,65],"A Helix",'g']]
    colors_regions = [regions_nm[0],regions_nm[-1]]                        
    # slice the regions 
    regions = [FEC_Util.slice_by_separation(example_plot,*reg) 
               for reg in zoom_regions_nm]
    ylim_pN_zoom = [50,120]
    # # make the plot 
    fig = PlotUtilities.figure((7,3))
    # create the 'top' gridspec
    top_spec = gridspec.GridSpec(1,2)
    # create separate axes for the image and FECs 
    image_spec = gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=top_spec[0])
    data_spec = gridspec.GridSpecFromSubplotSpec(3,1,subplot_spec=top_spec[1:],
                                                hspace=0.1) 
    # create a separate axis for the 'zoom in' FEC 
    zoom_in_spec = gridspec.GridSpecFromSubplotSpec(2,len(zoom_regions_nm),
                                                    subplot_spec=data_spec[1:],
                                                    wspace=0.06) 
    # # plot the image 
    ax = plt.subplot(image_spec[:])
    plt.imshow(plt.imread("../Data/sample_cartoon.png"),aspect='auto')
    ax.axis('off')
    xy_text = 0,np.mean(plt.ylim())
    # # plot the example fec and zoomed regions
    #
    # 'full' example 
    ax_example = plt.subplot(data_spec[0,:])
    alpha_data = 0.4
    color_data = 'k'
    dict_plot = dict(n_filter_points=2000,
                     style_data=dict(color=color_data,alpha=alpha_data,
                                     linewidth=0.5,linestyle='-'))
    x_full_plot = x_func(example_plot)                                     
    FEC_Plot._fec_base_plot(x_full_plot,y_func(example_plot),
                            **dict_plot)
    for i,(r,color) in enumerate(zip(regions,colors_regions)):
        dict_tmp = dict(**dict_plot)
        dict_tmp['style_data']['color'] = color[-1]
        FEC_Plot._fec_base_plot(x_func(r),y_func(r),**dict_tmp)
    plt.ylim(ylim_pN)
    plt.xlim(xlim_nm)
    ymin = 0.80
    ymax = 0.85
    for x,name,color in regions_nm:
        plt.axvspan(*x,color=color,ymin=ymin,ymax=ymax,alpha=0.3,linewidth=0)
    # plot the adhesion regions
    plt.axvspan(min(x_full_plot),adhesion_max_nm,color='k',alpha=0.1,
                hatch='//',linewidth=0)
    PlotUtilities.lazyLabel("Molecular extension (nm)","Force (pN)","")   
    PlotUtilities.x_label_on_top(ax_example)
    # # plot all the zoomed regions 
    offsets_x = [0.5,0.5]
    offsets_y = [0.7,0.15]
    heights_pN = [10,10]
    widths_nm = [1,1]
    for i,(r,color) in enumerate(zip(regions,colors_regions)):
        ax_tmp = plt.subplot(zoom_in_spec[:,i])
        dict_tmp = dict(**dict_plot)
        dict_tmp['style_data']['color'] = color[-1]
        FEC_Plot._fec_base_plot(x_func(r),y_func(r),**dict_tmp)
        xlim = plt.xlim()
        ylim = plt.ylim()
        min_x = min(xlim)
        max_y = max(ylim)
        y_loc = max_y*0.9
        kw = dict(ax=ax_tmp)
        x_kwargs =dict(unit="nm",width=widths_nm[i],**kw)
        offset_x = Scalebar.rel_to_abs(ax_tmp,offsets_x[i],True)
        offset_y = Scalebar.rel_to_abs(ax_tmp,offsets_y[i],False)
        Scalebar.crossed_x_and_y(offset_x=offset_x,
                                 offset_y=offset_y,
                                 x_kwargs=x_kwargs,
                                 y_kwargs=dict(unit="pN",height=heights_pN[i],
                                               **kw))
        PlotUtilities.no_x_label(ax_tmp)
        PlotUtilities.no_y_label(ax_tmp)        
        PlotUtilities.zoom_effect01(ax_example, ax_tmp,*xlim)

        PlotUtilities.lazyLabel("","","")           
    loc = [ [0.15,1.0],
            [-0.05,1.15],
            [0.1,0.9],
            [0.3,0.9]]
    PlotUtilities.label_tom(fig,loc=loc)
    subplots_adjust = dict(bottom=0.05,top=0.9,left=0)
    PlotUtilities.save_png_and_svg(fig,"diagram",
                                   subplots_adjust=subplots_adjust)

    
if __name__ == "__main__":
    run()
