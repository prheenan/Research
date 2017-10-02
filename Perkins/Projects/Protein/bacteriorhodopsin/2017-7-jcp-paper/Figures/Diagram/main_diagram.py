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
sys.path.append("../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import CheckpointUtilities,PlotUtilities,GenUtilities
from GeneralUtil.python.Plot import Scalebar,Annotations
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao
import figure_recreation
from figure_recreation import fig1d,fig4ab,fig4c

import jcp_fig_util

import copy 
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.patches

def velocity_annotate(ax,v,y=0.9,x=0.8,color='g'):
    s = r"$v$ = {:d} nm/s".format(int(np.round(v)))
    Annotations.relative_annotate(ax=ax,s=s,
                                  xy=(x,y),
                                  color=color,bbox=dict(color='w',pad=0))
                           
def run():
    """
    <Description>q

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
    example_plot.LowResData.force -= 7.1
    plot_examples = [example_plot]
    vel_m_per_s = example_plot.Velocity
    x_func = lambda y: y.Separation
    y_func = lambda y: y.Force 
    ylim_pN = [-20,155]
    xlim_nm = [5,75]
    zoom_regions_nm = [ [61.5,63.6]]
    adhesion_max_nm = 19
    region_labels = ["ED Helix","CB Helix","A Helix"]
    region_colors = jcp_fig_util.regions_and_colors()
    regions_nm = [[x,l,c] for l,(x,c) in zip(region_labels,region_colors)]
    colors_regions = [regions_nm[-1]]                        
    # slice the regions 
    regions = [FEC_Util.slice_by_separation(example_plot,*reg) 
               for reg in zoom_regions_nm]
    ylim_pN_zoom = [50,120]
    # # make the plot 
    fig = PlotUtilities.figure((7,4))
    # create the 'top' gridspec
    top_spec = gridspec.GridSpec(2,3)
    # create separate axes for the image and FECs 
    image_spec = \
        gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=top_spec[0,0])
    data_spec = \
        gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=top_spec[0,1:])
    # # plot the image 
    ax = plt.subplot(image_spec[:])
    plt.imshow(plt.imread("../Data/sample_cartoon.png"),aspect='auto')
    ax.axis('off')
    # # plot the example fec and zoomed regions
    #
    # 'full' example 
    ax_example = plt.subplot(data_spec[:,0])
    alpha_data = 0.4
    color_data = 'g'
    dict_plot = dict(n_filter_points=2000,
                     style_data=dict(color=color_data,alpha=alpha_data,
                                     linewidth=0.5,linestyle='-'))
    x_full_plot = x_func(example_plot)                                     
    FEC_Plot._fec_base_plot(x_full_plot,y_func(example_plot),
                            **dict_plot)
    PlotUtilities.tom_ticks(ax=ax_example,num_major=5,change_x=False)    
    PlotUtilities.tom_ticks(ax=ax_example,num_major=4,change_y=False)                              
    for i,(r,color) in enumerate(zip(regions,colors_regions)):
        # put a box around the region 
        x,y = x_func(r),y_func(r)
        Annotations.add_rectangle(ax_example,[min(x),max(x)],[min(y),max(y)])
    plt.ylim(ylim_pN)
    plt.xlim(xlim_nm)
    jcp_fig_util.add_helical_boxes(ax=ax_example)
    # plot the adhesion regions
    plt.axvspan(min(x_full_plot),adhesion_max_nm,color='0.85',
                linewidth=0)
    PlotUtilities.lazyLabel("Extension","Force","")   
    PlotUtilities.x_label_on_top(ax_example)
    PlotUtilities.no_x_label(ax_example)
    PlotUtilities.no_y_label(ax_example)
    x_kwargs = dict(unit_kwargs=dict(fmt="{:.0f}"),width=15,unit="nm")
    y_kwargs = dict(unit_kwargs=dict(fmt="{:.0f}"),
                    height=40,unit="pN")
    Scalebar.crossed_x_and_y_relative(offset_x=0.55,offset_y=0.58,
                                      x_kwargs=x_kwargs,
                                      y_kwargs=y_kwargs,
                                      ax=ax_example)    
    # add in the velocity annotation (in nm/s, from m/s)
    velocity_annotate(ax=ax_example,v=vel_m_per_s*1e9)
    # # plot all the zoomed regions 
    offsets_x = [0.8]
    offsets_y = [0.67]
    heights_pN = [10]
    widths_s = [0.001]
    for i,(r,color) in enumerate(zip(regions,colors_regions)):
        ax_tmp = plt.subplot(data_spec[-1])
        dict_tmp = dict(**dict_plot)
        dict_tmp['style_data']['color'] = color[-1]
        time = r.Time
        force = y_func(r)
        FEC_Plot._fec_base_plot(time,force,**dict_tmp)
        # add eyeballed lines to the data
        slope_1 = 10/3e-3
        slope_2 = 10/3.6e-3
        offset_1 = force[0] + 1 
        offset_2 = force[0] - 16.5
        line_upper = (time - min(time)) * slope_1 + offset_1
        line_lower = (time - min(time)) * slope_2 + offset_2
        time_rel = time - min(time)
        # dont let the line cross the textbox
        x_upper_idx = np.where(time_rel < 0.67 * max(time_rel))
        plt.plot(time[x_upper_idx],line_upper[x_upper_idx],'b--')
        plt.plot(time,line_lower,'b--')
        xlim = plt.xlim()
        ylim = plt.ylim()
        # add a little to the x and xlims, so the scale bar has a natural place
        # place to go 
        y_range = abs(np.diff(ylim)) 
        x_range = abs(np.diff(xlim)) 
        ylim = [ylim[0],ylim[1]+y_range*0.05]
        xlim = [xlim[0],xlim[1]+x_range*0.01]        
        ax_tmp.set_ylim(ylim)
        ax_tmp.set_xlim(xlim)
        min_x = min(xlim)
        max_y = max(ylim)
        y_loc = max_y*0.9
        x_kwargs =dict(unit="ms",width=widths_s[i],
                       fudge_text_pct=dict(x=0.2,y=0),
                       unit_kwargs=dict(value_function=lambda x: x * 1e3))
        y_kwargs = dict(unit="pN ",
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
    # read in the data ...
    base_re = "./recreation_figure_data/"
    # # Figure 1D from the science paper 
    ax_fec_ensemble = plt.subplot(top_spec[1,0])
    fig4d = figure_recreation.save_output(base_re,"Fig1D.csv")      
    ylim = [0,160]
    xlim = [18,32]
    for wlc_x,wlc_y in zip(fig4d.wlc_x,fig4d.wlc_y):
        plt.plot(wlc_x,wlc_y,'b--',alpha=0.4,linewidth=1,dashes=(2,2))
    for x,y in zip(fig4d.x,fig4d.y):
        plt.plot(x,y,alpha=1,linewidth=0.5)        
    ax_fec_ensemble.set_ylim(ylim)
    ax_fec_ensemble.set_xlim(xlim)    
    PlotUtilities.lazyLabel("Extension","Force","")        
    x_kwargs = dict(width=3,unit="nm")
    y_font = copy.deepcopy(Scalebar.def_font_kwargs_y)
    y_font['rotation'] = 90
    y_kwargs = dict(height=25,unit="pN",font_kwargs=y_font)
    Scalebar.crossed_x_and_y_relative(offset_x=0.22,offset_y=0.77,
                                      x_kwargs=x_kwargs,
                                      y_kwargs=y_kwargs,
                                      ax=ax_fec_ensemble)    
    PlotUtilities.no_x_label(ax=ax_fec_ensemble)                                      
    PlotUtilities.no_y_label(ax=ax_fec_ensemble)  
    # # Figure 4B from the science paper
    color_equil = 'rebeccapurple'
    fig4ab = figure_recreation.save_output(base_re,"Fig4AB.csv")  
    ax_time = plt.subplot(top_spec[1,1])
    min_x,max_x = min(fig4ab.time),max(fig4ab.time)
    range = [0.55,0.59]
    min_x_new = min_x + (max_x-min_x)*range[0]
    max_x_new = min_x + (max_x-min_x)*range[1]
    idx = np.where( (fig4ab.time <= max_x_new) & (fig4ab.time >= min_x_new))
    time = fig4ab.time[idx]
    force = fig4ab.force[idx]
    FEC_Plot._fec_base_plot(time,force,n_filter_points=200,
                            style_data=dict(color=color_equil,alpha=0.3,
                                            linewidth=0.75))
    ax_time.set_xlim(min_x_new,max_x_new)
    ax_time.set_ylim(None,None)
    unit_kwargs = dict(value_function =lambda x: x*1e6,fmt="{:.0f}")
    unit_micro_s = PlotUtilities.upright_mu() + "m"
    x_kwargs = dict(unit_kwargs=unit_kwargs,width=500e-6,unit=unit_micro_s)
    y_font = copy.deepcopy(Scalebar.def_font_kwargs_y)
    y_font['rotation'] = 90
    y_kwargs = dict(height=10,unit="pN",font_kwargs=y_font)
    Scalebar.crossed_x_and_y_relative(offset_x=0.5,offset_y=0.08,
                                      x_kwargs=x_kwargs,
                                      y_kwargs=y_kwargs,
                                      ax=ax_time)  
    PlotUtilities.no_x_label(ax=ax_time)                                      
    PlotUtilities.no_y_label(ax=ax_time)  
    PlotUtilities.lazyLabel("Time","Force","")
    velocity_annotate(ax=ax_time,v=0,color=color_equil)
    # # figure 4C from the science paper -- the pfold energy landscape 
    fig4c = figure_recreation.save_output(base_re,"Fig4C.csv")
    ax_equil = plt.subplot(top_spec[1,2])
    # data is in kJ/mol, communication with hao, 2017-9-14
    fig4c.energy /= 4.2
    fig4c.energy_error /= 4.2
    plt.errorbar(fig4c.x,fig4c.energy,fig4c.energy_error,color=color_equil,
                 marker='o',
                 mfc='w',zorder=0,markerfacecolor="None",capsize=2,elinewidth=1,
                 linewidth=1)                 
    PlotUtilities.lazyLabel("Extension (nm)","Energy","")     
    x_kwargs = dict(unit_kwargs=dict(fmt="{:.1f}"),width=0.1,unit="nm")
    y_font = copy.deepcopy(Scalebar.def_font_kwargs_y)
    y_font['rotation'] = 90
    y_kwargs = dict(unit_kwargs=dict(fmt="{:.1f}"),
                    height=0.5,unit="kcal/mol",font_kwargs=y_font)
    Scalebar.crossed_x_and_y_relative(offset_x=0.22,offset_y=0.60,
                                      x_kwargs=x_kwargs,
                                      y_kwargs=y_kwargs,
                                      ax=ax_equil)
    # add in bell...
    bell_mean = 0.64
    bell_std = 0.09
    x0,xf = ax_equil.get_xlim()
    mean_x = bell_mean * (xf-x0) + x0
    std_x = bell_std * (xf-x0)
    color_bell = 'k'
    plt.axvspan(mean_x-std_x,mean_x+std_x,color=color_bell,alpha=0.15)
    plt.axvline(mean_x,color=color_bell,linestyle='--',zorder=0,alpha=0.7)
    PlotUtilities.no_x_label(ax=ax_equil)         
    t = ax_equil.annotate(s=r"$\Delta x^{\ddag}_{\mathrm{Bell}}$",
                          xy=(0.4,0.3),color=color_bell,
                          xycoords="axes fraction")
    PlotUtilities.no_y_label(ax=ax_equil)    
    loc_upper = [-0.05,1.05]
    loc_lower = [-0.05,1.0]
    loc = [loc_upper,loc_upper,loc_upper,   
           loc_lower,loc_lower,loc_lower]
    PlotUtilities.label_tom(fig,loc=loc)
    subplots_adjust = dict(hspace=0.07,wspace=0.15)
    PlotUtilities.save_png_and_svg(fig,"diagram",
                                   subplots_adjust=subplots_adjust)
    
    
if __name__ == "__main__":
    run()
