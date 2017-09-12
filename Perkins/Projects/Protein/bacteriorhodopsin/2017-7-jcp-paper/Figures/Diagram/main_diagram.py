# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,matplotlib as mpl


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
import figure_recreation
from figure_recreation import fig1d,fig4ab,fig4c


import copy 
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.patches

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = None
    flickering_dir = "../Data/fec/"
    # XXX use the flickering dir for stuff
    cache_dir = flickering_dir 
    GenUtilities.ensureDirExists(flickering_dir)
    force_read_data = False
    raw_data = IoUtilHao.read_and_cache_data_hao(flickering_dir,
                                                force=force_read_data,
                                                 cache_directory=flickering_dir,
                                                 limit=3)
    example = raw_data[0]                                                 
    example_plot = copy.deepcopy(example)
    # fix the manual offset
    example_plot.Force -= 7.1
    plot_examples = [example_plot]
    x_func = lambda y: y.Separation
    y_func = lambda y: y.Force 
    ylim_pN = [-20,None]
    xlim_nm = [-5,100]
    zoom_regions_nm = [ [61.5,63]]
    adhesion_max_nm = 19
    # plot the helical regions...
    regions_nm = [ [[adhesion_max_nm,30],"ED Helix",'royalblue'],
                   [[31,45],"CB Helix",'orangered'],
                   [[55,65],"A Helix",'g']]
    colors_regions = [regions_nm[-1]]                        
    # slice the regions 
    regions = [FEC_Util.slice_by_separation(example_plot,*reg) 
               for reg in zoom_regions_nm]
    ylim_pN_zoom = [50,120]
    # # make the plot 
    fig = PlotUtilities.figure((7,4))
    # create the 'top' gridspec
    top_spec = gridspec.GridSpec(2,3,left=0)
    # create separate axes for the image and FECs 
    image_spec = \
        gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=top_spec[0,0],
                                        hspace=0,
                                        wspace=0)
    data_spec = \
        gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=top_spec[0,1:],
                                         hspace=0.05,wspace=0.12)
    # # plot the image 
    ax = plt.subplot(image_spec[:])
    plt.imshow(plt.imread("../Data/sample_cartoon.png"),aspect='auto')
    ax.axis('off')
    # # plot the example fec and zoomed regions
    #
    # 'full' example 
    ax_example = plt.subplot(data_spec[:,0])
    alpha_data = 0.4
    color_data = 'k'
    dict_plot = dict(n_filter_points=1000,
                     style_data=dict(color=color_data,alpha=alpha_data,
                                     linewidth=0.5,linestyle='-'))
    x_full_plot = x_func(example_plot)                                     
    FEC_Plot._fec_base_plot(x_full_plot,y_func(example_plot),
                            **dict_plot)
    PlotUtilities.tom_ticks(ax=ax_example,num_major=5,change_x=False)    
    PlotUtilities.tom_ticks(ax=ax_example,num_major=4,change_y=False)                              

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
    mpl.rcParams['hatch.color'] = '0.85'
    plt.axvspan(min(x_full_plot),adhesion_max_nm,color='0.95',hatch='--',
                linewidth=0)
    PlotUtilities.lazyLabel("Molecular extension (nm)","Force (pN)","")   
    PlotUtilities.x_label_on_top(ax_example)
    # # plot all the zoomed regions 
    offsets_x = [0.9,0.5]
    offsets_y = [0.7,0.15]
    heights_pN = [5,10]
    widths_s = [0.001]
    for i,(r,color) in enumerate(zip(regions,colors_regions)):
        ax_tmp = plt.subplot(data_spec[-1])
        dict_tmp = dict(**dict_plot)
        dict_tmp['style_data']['color'] = color[-1]
        FEC_Plot._fec_base_plot(r.Time,y_func(r),**dict_tmp)
        xlim = plt.xlim()
        ylim = plt.ylim()
        min_x = min(xlim)
        max_y = max(ylim)
        y_loc = max_y*0.9
        line_kw = dict(linewidth=1.0,color='k')
        x_kwargs =dict(unit="ms",width=widths_s[i],line_kwargs=line_kw,
                       unit_kwargs=dict(value_function=lambda x: x * 1e3))
        y_kwargs = dict(unit="pN ",line_kwargs=line_kw,
                        height=heights_pN[i])
        offset_x = Scalebar.rel_to_abs(ax_tmp,offsets_x[i],True)
        offset_y = Scalebar.rel_to_abs(ax_tmp,offsets_y[i],False)
        Scalebar.crossed_x_and_y(offset_x=offset_x,
                                 offset_y=offset_y,
                                 ax=ax_tmp,
                                 x_kwargs=x_kwargs,
                                 y_kwargs=y_kwargs)
        PlotUtilities.no_y_label(ax_tmp)
        PlotUtilities.no_x_label(ax_tmp)
        PlotUtilities.x_label_on_top(ax_tmp)        
        PlotUtilities.lazyLabel("","","")        
        PlotUtilities.xlabel("Time")
    # add a single arrow on the last axis. See:
    # https://www.cilyan.org/blog/2016/01/23/matplotlib-draw-between-subplots/
    # 1. Get transformation operators for axis and figure
    ax0tr = ax_example.transData # Axis 0 -> Display
    ax1tr = ax_tmp.transData # Axis 1 -> Display
    figtr = fig.transFigure.inverted() # Display -> Figure
    # 2. Transform arrow start point from axis 0 to figure coordinates
    ptB = figtr.transform(ax0tr.transform((85,50)))
    # 3. Transform arrow end point from axis 1 to figure coordinates
    tmp_ylim = ax_tmp.get_ylim()
    max_y = np.max(tmp_ylim)
    range_y = abs(np.diff(tmp_ylim))
    tmp_xlim = ax_tmp.get_xlim()
    range_x = abs(np.diff(tmp_xlim))
    y_styles = [ (max_y-range_y*0.9,"arc3,rad=0.2"),
                 (max_y-range_y*0.3,"arc3,rad=-0.2")]
    for y,style in y_styles:                  
        ptE = figtr.transform(ax1tr.transform((np.mean(tmp_xlim)*0.999,
                                               y)))                                           
        # 4. Create the patch
        arrow = matplotlib.patches.FancyArrowPatch(
            ptB, ptE, transform=fig.transFigure,  # Place arrow in figure coord system
            fc = "g", connectionstyle=style, arrowstyle='simple', 
            alpha = 0.3,linewidth=0,
            mutation_scale = 20.
        )
        # 5. Add patch to list of objects to draw onto the figure
        fig.patches.append(arrow)
    # read in the data ...
    base_re = "./recreation_figure_data/"
    # # Figure 4B from the science paper
    fig4ab = figure_recreation.save_output(base_re,"Fig4AB.csv")  
    ax_time = plt.subplot(top_spec[1,1])
    plt.plot(fig4ab.time,fig4ab.force,linewidth=0.5)
    min_x,max_x = plt.xlim()
    range = [0.55,0.61]
    min_x_new = min_x + (max_x-min_x)*range[0]
    max_x_new = min_x + (max_x-min_x)*range[1]
    ax_time.set_xlim(min_x_new,max_x_new)
    ax_time.set_ylim(None,None)
    unit_kwargs = dict(value_function =lambda x: x*1000,fmt="{:.0f}")
    x_kwargs = dict(unit_kwargs=unit_kwargs,width=0.001,unit="ms")
    y_font = copy.deepcopy(Scalebar.def_font_kwargs_y)
    y_font['rotation'] = 90
    y_kwargs = dict(height=10,unit="pN",font_kwargs=y_font)
    Scalebar.crossed_x_and_y_relative(offset_x=0.7,offset_y=0.08,
                                      x_kwargs=x_kwargs,
                                      y_kwargs=y_kwargs,
                                      ax=ax_time)  
    PlotUtilities.no_x_label(ax=ax_time)                                      
    PlotUtilities.no_y_label(ax=ax_time)  
    PlotUtilities.lazyLabel("Time","Force","")
    # # figure 4C from the science paper -- the pfold energy landscape 
    fig4c = figure_recreation.save_output(base_re,"Fig4C.csv")
    ax_equil = plt.subplot(top_spec[1,2])
    plt.errorbar(fig4c.x,fig4c.energy,fig4c.energy_error,fmt='ko-',
                 mfc='w',zorder=0)                 
    PlotUtilities.lazyLabel("Extension (nm)","Energy (kcal/mol)","")     
    x_kwargs = dict(unit_kwargs=dict(fmt="{:.1f}"),width=0.1,unit="nm")
    y_font = copy.deepcopy(Scalebar.def_font_kwargs_y)
    y_font['rotation'] = 90
    y_kwargs = dict(height=2,unit="kcal/mol",font_kwargs=y_font)
    Scalebar.crossed_x_and_y_relative(offset_x=0.25,offset_y=0.65,
                                      x_kwargs=x_kwargs,
                                      y_kwargs=y_kwargs,
                                      ax=ax_equil)
    PlotUtilities.no_x_label(ax=ax_equil)     
    # XXX figure out what is wrong with this?
    #PlotUtilities.no_y_label(ax=ax_equil)                                          
    loc = [ [0.15,1.15],
            [-0.05,1.15],
            [-0.05,1.15]]
    PlotUtilities.label_tom(fig,loc=loc)
    PlotUtilities.save_png_and_svg(fig,"diagram")
    
if __name__ == "__main__":
    run()
