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

sys.path.append("../")
sys.path.append("../../../../../../../../")
import jcp_fig_util,GenerateLandscapes
from GenerateLandscapes import slice_area,cacheable_data
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import CheckpointUtilities,PlotUtilities,GenUtilities
from GeneralUtil.python.Plot import Scalebar,Annotations 
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass,WeierstrassUtil
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


import copy 
from matplotlib import gridspec
import scipy

def make_heatmap(histogram, x_edges,y_edges,kw_heatmap):
    # XXX ? digitize all the ids so we know what bin they fall into...
    X,Y = np.meshgrid(x_edges,y_edges)
    X -= np.min(X)
    plt.gca().pcolormesh(X,Y,histogram,**kw_heatmap)
        
def plot_landscape(data,xlim,kw_landscape=dict(),plot_derivative=True,
                   zero_q=True,
                   label_deltaG = PlotUtilities.variable_string("\Delta G")):
    landscape_kcal_per_mol = data.mean_landscape_kcal_per_mol
    landscape_kcal_per_mol -= min(landscape_kcal_per_mol)
    std_landscape_kcal_per_mol = data.std_landscape_kcal_per_mol
    extension_nm = data._extension_grid_nm
    if (zero_q):
        extension_nm -= min(extension_nm)
    extension_aa = data.amino_acids_per_nm() * extension_nm
    grad = lambda x: np.gradient(x)/(np.gradient(extension_aa))
    delta_landscape_kcal_per_mol_per_amino_acid = grad(landscape_kcal_per_mol)
    landscape_upper = landscape_kcal_per_mol+std_landscape_kcal_per_mol
    landscape_lower =landscape_kcal_per_mol-std_landscape_kcal_per_mol
    # XXX -- everything is in nm now, by tom's request. should refator
    delta_landscape_kcal_per_mol_per_amino_acid = \
        data.mean_delta_landscape_kcal_per_mol_per_nm
    std_delta_landscape_kcal_per_mol_per_AA = \
        data.std_delta_landscape_kcal_per_mol_per_nm
    upper_delta_landscape = delta_landscape_kcal_per_mol_per_amino_acid+\
                            std_delta_landscape_kcal_per_mol_per_AA
    lower_delta_landscape = delta_landscape_kcal_per_mol_per_amino_acid-\
                            std_delta_landscape_kcal_per_mol_per_AA
    # make a second axis for the number of ammino acids 
    limits_delta = [min(delta_landscape_kcal_per_mol_per_amino_acid),
                    max(delta_landscape_kcal_per_mol_per_amino_acid)]    
    limits_energy = [min(landscape_kcal_per_mol),max(landscape_kcal_per_mol)]
    units_energy_delta = label_deltaG + \
                         r" per nm (kcal/(mol $\cdot$ nm))"   
    units_energy = PlotUtilities.unit_string("\Delta G","kcal/mol")                         
    difference_color = 'rebeccapurple'    
    landscape_color = kw_landscape['color']
    if (plot_derivative):
        PlotUtilities.lazyLabel("Extension (nm)",units_energy_delta,"")
        PlotUtilities.ylabel(units_energy_delta,color=difference_color)
    else:
        PlotUtilities.lazyLabel("Extension (nm)","","")  
        PlotUtilities.ylabel(units_energy,color=landscape_color)       
    # the derivative goes on the main axis in this case...
    if (plot_derivative):
        ax_delta = plt.gca()
        ax_2 = PlotUtilities.secondAxis(ax_delta,
                                        label=units_energy,
                                        color=landscape_color,
                                        limits=limits_energy,secondY =True,
                                        tick_color=landscape_color)
        ax_energy = ax_2
        to_ret = ax_delta,ax_energy
    else:
        ax_energy = plt.gca()
        ax_delta = None
        to_ret = ax_energy
    # filter the data if we need to 
    if (xlim is not None and xlim[1] is not None):
        idx= np.where( (extension_nm > xlim[0]) & (extension_nm < xlim[1]))
    else:
        idx =np.arange(extension_nm.size)
    # plot the landscape and its standard deviation
    ax_energy.plot(extension_nm[idx],landscape_kcal_per_mol[idx],
                   label=label_deltaG,zorder=10,**kw_landscape)
    ax_energy.fill_between(x=extension_nm[idx],
                           y1=landscape_lower[idx],
                           y2=landscape_upper[idx],
                           alpha=0.15,zorder=10,**kw_landscape)      
    # the energy y label should be rotated if it is on the right                            
    rotation = -90    
    if plot_derivative:
        kw_y = dict(rotation=-90,x=-0.3)
    else:
        kw_y = dict()
    PlotUtilities.ylabel(ax=ax_energy,lab=units_energy,**kw_y)
    # move the y label to the right slightly i we are using both axes 
    if plot_derivative:
        ax_energy.yaxis.set_label_coords(1.2,0.5)
    if (plot_derivative):                           
        # plot the energy delta and its bounds, based on the bounds on the
        #        landscape    
        ax_delta.plot(extension_nm[idx],
                      delta_landscape_kcal_per_mol_per_amino_acid[idx],
                      color=difference_color,linestyle='-',linewidth=1.5)
        PlotUtilities.color_axis_ticks(color=landscape_color,ax=ax_energy,
                                       spine_name='right')   
        ax_delta.fill_between(x=extension_nm[idx],
                              y1=lower_delta_landscape[idx],
                              y2=upper_delta_landscape[idx],
                              alpha=0.15,color=difference_color)
        PlotUtilities.tom_ticks(ax=ax_delta,change_x=False,num_major=4)
    return to_ret
              
def heatmap_plot(heatmap_data,amino_acids_per_nm,kw_heatmap=dict()):
    ax_heat = plt.gca()
    make_heatmap(*heatmap_data,kw_heatmap=kw_heatmap)
    PlotUtilities.lazyLabel("","Force (pN)","")
    # make a second x axis for the number of ammino acids 
    xlim_fec = plt.xlim()
    limits = np.array(xlim_fec) * amino_acids_per_nm
    tick_kwargs = dict(axis='both',color='w',which='both')                                 
    ax_heat.tick_params(**tick_kwargs)                                         
    plt.xlim(xlim_fec)
    
def create_landscape_plot(data_to_plot,kw_heatmap=dict(),kw_landscape=dict(),
                          xlim=None,zero_q=True):
    """
    Creates a plot of
    """
    heatmap_data = data_to_plot.heatmap_data
    data_landscape = data_to_plot.generate_landscape_obj()
    # # ploy the heat map 
    ax_heat = plt.subplot(2,1,1)
    heatmap_plot(heatmap_data,data_landscape.amino_acids_per_nm(),
                 kw_heatmap=kw_heatmap)
    # # plot the energy landscape...
    ax_energy = plt.subplot(2,1,2)    
    ax1,ax2 = plot_landscape(data_landscape,xlim,kw_landscape=kw_landscape,
                             zero_q=zero_q)
    if (xlim is None):
        xlim = ax1.get_xlim()
        xlim = np.maximum(0,xlim)
    normalize_and_set_zeros([ax1],[ax2])
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)
    ax_heat.set_xlim(xlim)
    return ax_heat,ax_energy

def landscape_kwargs():
    """
    Returns: a list of keywords for the entire listscape, ED,CB, and A helices 
    """
    kwargs = [ dict(kw_heatmap=dict(cmap=plt.cm.afmhot),
                    kw_landscape=dict(color='k')),
               dict(kw_heatmap=dict(cmap=plt.cm.Blues_r),
                    kw_landscape=dict(color='royalblue')),
               dict(kw_heatmap=dict(cmap=plt.cm.Reds_r),
                    kw_landscape=dict(color='orangered')),                    
               dict(kw_heatmap=dict(cmap=plt.cm.Greens_r),
                    kw_landscape=dict(color='g'))]
    return kwargs                    
    
def make_detalied_plots(data_to_analyze,areas):
    """
    makes the detailed plots
    """
    kwargs = landscape_kwargs()
    y_limits_pN = [None,None,125,100]
    for i,a in enumerate(areas):
        fig = PlotUtilities.figure((3.25,5))
        mdata = data_to_analyze[i]
        example = mdata.landscape[0]
        ax = create_landscape_plot(mdata,xlim=None,zero_q=False,**kwargs[i])
        ax_heat = ax[0]
        PlotUtilities.no_x_label(ax_heat)
        ax_heat.relim()
        ax_heat.set_ylim([None,y_limits_pN[i]])
        out_name = "landscape{:d}_{:s}".format(i,areas[i].plot_title)
        axis_func = lambda x: [x[0],x[2]]
        PlotUtilities.label_tom(fig,axis_func=axis_func)
        PlotUtilities.save_png_and_svg(fig,out_name.replace(" ","_"))    
    
def normalize_axes(ax_list,manual_min=None,fudge_f=0.0):
    """
    ensures all the y axes in ax_list have the same limits
    """
    max_ylim = np.max([ax.get_ylim() for ax in ax_list])
    min_ylim = np.min([ax.get_ylim() for ax in ax_list])
    range_v = max_ylim-min_ylim
    fudge = range_v * fudge_f
    min_v = min_ylim-fudge
    if (manual_min is not None):
        min_v = manual_min    
    for ax in ax_list:
        ax.set_ylim([min_v,max_ylim+fudge/2])    
    
def helical_gallery_plot(helical_areas,helical_data,helical_kwargs):
    axs,first_axs,second_axs = [],[],[]
    offset_y = 0.2
    offset_x = 0.30
    kw_scalebars = [dict(offset_x=offset_x,offset_y=offset_y),
                    dict(offset_x=offset_x,offset_y=offset_y),
                    dict(offset_x=offset_x,offset_y=offset_y)]
    xlims = [ [None,None],[None,None],[None,15]    ]   
    arrow_x = [0.65,0.67,0.47]
    arrow_y = [0.60,0.52,0.42]
    for i,a in enumerate(helical_areas):
        data = helical_data[i]
        kw_tmp = helical_kwargs[i]
        data_landscape = data.generate_landscape_obj()
        # # plot the energy landscape...
        ax_tmp = plt.subplot(1,len(helical_areas),(i+1))
        axs.append(ax_tmp)
        kw_landscape = kw_tmp['kw_landscape']
        color = kw_landscape['color']      
        ax_1, ax_2 = plot_landscape(data_landscape,xlim=xlims[i],
                                    kw_landscape=kw_landscape,
                                    plot_derivative=True)
        first_axs.append(ax_1)                                 
        second_axs.append(ax_2)                    
        PlotUtilities.tom_ticks(ax=ax_2,num_major=5,change_x=False)       
        last_idx = len(helical_areas)-1
        ax_1.annotate("",xytext=(arrow_x[i],arrow_y[i]),
                      textcoords='axes fraction',
                      xy=(arrow_x[i]+0.2,arrow_y[i]),xycoords='axes fraction',
                      arrowprops=dict(facecolor=color,alpha=0.7,
                                      edgecolor="None",width=4,headwidth=10,
                                      headlength=5))
        if (i > 0):
            PlotUtilities.ylabel("")
            PlotUtilities.xlabel("")
        if (i != last_idx):
            ax_2.set_ylabel("")
            PlotUtilities.no_x_label(ax_1)
            PlotUtilities.no_x_label(ax_2)            
        PlotUtilities.title(a.plot_title,color=color,y=0.97)
    normalize_and_set_zeros(first_axs,second_axs)
    # after normalization, add in the scale bars 
    for i,(ax_1,ax_2) in enumerate(zip(first_axs,second_axs)):
        Scalebar.x_scale_bar_and_ticks_relative(unit="nm",width=5,ax=ax_2,
                                                **kw_scalebars[i])
        PlotUtilities.no_x_label(ax_2)     
            
    
def normalize_and_set_zeros(first_axs,second_axs,fudge_1=0.05,fudge_2=0.6):
    # make all the axes consistent            
    normalize_axes(first_axs,fudge_f=fudge_1)
    normalize_axes(second_axs,fudge_f=fudge_2)
    # determine where the min need to be for the zeros to be OK
    min_y,max_y = first_axs[0].get_ylim()
    normalized_zero = -min_y/(max_y-min_y)
    min_y_second,max_y_second = second_axs[0].get_ylim()
    # want f = 0-min/(max-min) --> min = f*max/(f-1)
    needed_min = normalized_zero * max_y_second / (normalized_zero-1)
    normalize_axes(second_axs,manual_min=needed_min)       
    
def make_gallery_plot(areas,data_to_analyze,out_name="./gallery"):
    # skip the first one (the entire landscape )
    helical_areas = areas[1:]
    helical_data = data_to_analyze[1:]
    helical_kwargs = landscape_kwargs()[1:]
    fig = PlotUtilities.figure((7,2.5))
    helical_gallery_plot(helical_areas,helical_data,helical_kwargs)    
    ax_to_label = [ax for i,ax in enumerate(fig.axes) if (i % 2) != 0]
    loc_last_two = [-0.1,1.05]
    loc_first = [-0.2,1.05]
    PlotUtilities.label_tom(fig=fig,axis_func=lambda *a,**kw: ax_to_label,
                            loc=[loc_first,loc_last_two,loc_last_two])
    PlotUtilities.save_png_and_svg(fig,out_name)                    
    
def setup_pedagogy_ticks(ax,scale_bar_x,x_heat_kw,y_heat_kw,offset_y=0.9):
    font_kwargs= copy.deepcopy(Scalebar.def_font_kwargs_y)
    # fix the keywords relative to the heatmap
    x_heat_kw['font_kwargs']['color'] = 'k'
    x_heat_kw['line_kwargs']['color'] = 'k'
    y_heat_kw['font_kwargs']['color'] = 'k'
    y_heat_kw['line_kwargs']['color'] = 'k'  
    y_heat_kw['unit'] = 'kcal/mol '
    y_heat_kw['height'] = 50
    Scalebar.crossed_x_and_y_relative(scale_bar_x,offset_y,ax=ax,
                                      x_kwargs=x_heat_kw,
                                      y_kwargs=y_heat_kw)
    PlotUtilities.no_y_label(ax)            
    # add in a zero point
    fudge_x = 1
    xlim = ax.get_xlim()
    fontsize = 7
    ax.annotate(s="0",
                xytext=(xlim[0] - fudge_x,0),
                xy=(0,0),
                size=fontsize,
                textcoords='data',
                verticalalignment="center",
                horizontalalignment="center")


def kwargs_correction():     
    return [dict(color='teal',linestyle='-.'),
            dict(color='m',linestyle='-',linewidth=0.75),
            dict(color='saddlebrown',linestyle='--')]
            
def kwargs_labels():
    second_deriv =  r"\frac{1}{2\beta}\ln(1-\frac{\ddot{A}}{k})"                
    return [PlotUtilities.variable_string(r"A"),
            PlotUtilities.variable_string(r"\frac{\dot{A}^2}{2k}"),
            PlotUtilities.variable_string(second_deriv)]


            
def plot_with_corrections(data):
    ext_nm = 1e9 * data._extensions_m_original.copy()
    ext_nm -= min(ext_nm)
    convert = data.from_Joules_to_kcal_per_mol()
    key = data.landscape_objs[0]
    # the second derivative is a little special; we average the argument to
    # the log first,
    # to avoid lots of zeros
    beta = key.beta
    f_deriv = data._grid_property(lambda x: x.one_minus_A_z_ddot_over_k)
    f_mean = lambda x: np.mean(x,axis=0)
    to_second_deriv = lambda x: InverseWeierstrass.second_deriv_term(x,beta)
    first_deriv = \
        data._grid_property(lambda x: -1 * x.first_deriv_term * convert)
    energies = [f_mean(data._grid_property(lambda x: x.A_z * convert)),
                f_mean(first_deriv),
                to_second_deriv(f_mean(f_deriv))*convert]
    labels = kwargs_labels()
    kwargs = kwargs_correction()
    # plot each
    for i,energy_rel in enumerate(energies):
        plt.plot(ext_nm,energy_rel-min(energy_rel),label=labels[i],**kwargs[i])
    
def _second_deriv_plot(ax_heat,data):
    units_y = lambda x: x  * (1e12/1e9)
    mean_second_deriv_pN_nm = units_y(data._avg(data._d2_energies_dm2))
    stdev_second_deriv_pN_nm = units_y(data._std(data._d2_energies_dm2))
    q_nm = data._extension_grid_nm    
    upper_mean = mean_second_deriv_pN_nm + stdev_second_deriv_pN_nm
    lower_mean = mean_second_deriv_pN_nm - stdev_second_deriv_pN_nm
    # only plot the first xn; the A helix does something weird after detachment
    max_nm = 45
    where_to_plot = np.where(q_nm <= max_nm)                                   
    q_nm = q_nm[where_to_plot]
    upper_mean = upper_mean[where_to_plot]
    lower_mean = lower_mean[where_to_plot]
    mean_second_deriv_pN_nm = mean_second_deriv_pN_nm[where_to_plot]
    stdev_second_deriv_pN_nm = stdev_second_deriv_pN_nm[where_to_plot]
    cantilever_stiffness_pN_per_nm = 25
    label_cantlever = (r"$k_{\mathrm{cantilever}}\approx$" + \
                        "{:d} pN/nm".format(cantilever_stiffness_pN_per_nm))
    plt.axhline(cantilever_stiffness_pN_per_nm,color='m',linestyle='--',
                label=label_cantlever,linewidth=1.5)
    regions_and_colors = jcp_fig_util.regions_and_colors(subtract_min=True)
    for r,c in regions_and_colors:
        idx_tmp = np.where( (q_nm >= r[0]) & (q_nm <= r[1]))
        plt.plot(q_nm[idx_tmp],mean_second_deriv_pN_nm[idx_tmp],color=c,
                 linewidth=0.75)
        ax_heat.fill_between(x=q_nm[idx_tmp],
                             y1=lower_mean[idx_tmp],
                             y2=upper_mean[idx_tmp],
                             alpha=0.3,color=c,linewidth=0)  
    jcp_fig_util.add_helical_boxes(ax=ax_heat,ymax_box=0.1,alpha=1.0,
                                   font_color='w',offset_bool=True,
                                   max_x=max_nm) 
    PlotUtilities.lazyLabel("Extension (nm)","Stiffness (pN/nm)","",
                            legend_kwargs=dict(handlelength=2),loc='upper left')   
    
def make_second_deriv_plot(data_to_plot,kw):
    data = data_to_plot.generate_landscape_obj()
    # get the average and stdev of the second derivative
    # y is in units of N/m, convert to pN/nm (multiply by 1e3)

    fig = PlotUtilities.figure()
    ax_heat = plt.subplot(1,1,1)
    _second_deriv_plot(ax_heat,data)
    PlotUtilities.savefig(fig,"./out_stiff.png")
    
def make_pedagogical_plot(data_to_plot,kw,out_name="./iwt_diagram"):
    heatmap_data = data_to_plot.heatmap_data
    data = data_to_plot.generate_landscape_obj()
    fig = PlotUtilities.figure((3.25,5))
    # # ploy the heat map 
    ax_heat = plt.subplot(3,1,1)
    heatmap_plot(heatmap_data,data.amino_acids_per_nm(),
                 kw_heatmap=kw['kw_heatmap'])
    xlim_fec = plt.xlim()
    PlotUtilities.no_x_label(ax_heat)    
    ax_heat.set_ylim([0,150])
    PlotUtilities.no_x_label(ax_heat)
    PlotUtilities.no_y_label(ax_heat)   
    fontsize_scalebar = 6 
    common_kw = dict(color='w',fontsize=fontsize_scalebar)
    x_font,y_font = Scalebar.\
        font_kwargs_modified(x_kwargs=common_kw,
                             y_kwargs=common_kw)
    heat_kw_common = dict(line_kwargs=dict(color='w',linewidth=1.5))
    fudge_x = dict(x=0,y=0.25)
    fudge_y = dict(x=0,y=0)
    x_heat_kw = dict(width=15,unit="nm",font_kwargs=x_font,
                     fudge_text_pct=fudge_x,**heat_kw_common)
    y_heat_kw = dict(height=30,unit='pN ',
                     fudge_text_pct=fudge_y,font_kwargs=y_font,**heat_kw_common)
    # add a scale bar for the heatmap...
    scale_bar_x = 0.83
    Scalebar.crossed_x_and_y_relative(scale_bar_x,0.55,ax=ax_heat,
                                      x_kwargs=x_heat_kw,
                                      y_kwargs=y_heat_kw)
    jcp_fig_util.add_helical_boxes(ax=ax_heat,ymax_box=0.9,alpha=1.0,
                                   font_color='w',offset_bool=True)
    # # plot the energy landscape...
    ax_correction = plt.subplot(3,1,2)    
    plot_with_corrections(data)
    PlotUtilities.no_x_label(ax_correction)
    PlotUtilities.lazyLabel("","Energy (kcal/mol)","")
    ax_correction.set_xlim(xlim_fec)            
    offset_y_pedagogy = 0.31
    setup_pedagogy_ticks(ax_correction,scale_bar_x,x_heat_kw,y_heat_kw,
                         offset_y=offset_y_pedagogy)
    legend_font_size = 9                         
    legend = PlotUtilities.legend(handlelength=1.5,loc=(0.15,0.07),ncol=3,
                                  fontsize=legend_font_size,handletextpad=0.4)
    for i,text in enumerate(legend.get_texts()):
        plt.setp(text, color = kwargs_correction()[i]['color'])    
    # make the inset plot 
    axins = zoomed_inset_axes(ax_correction, zoom=2.25, loc=2,
                              borderpad=1.2) 
    plot_with_corrections(data)
    xlim_box = [-1,6]
    ylim_box = [-2,29]
    plt.xlim(xlim_box)
    plt.ylim(ylim_box)
    PlotUtilities.no_x_anything(axins)
    PlotUtilities.no_y_anything(axins)
    # add in scale bars
    kw_common = dict(line_kwargs=dict(linewidth=0.75,color='k'))
    common_font_inset = dict(fontsize=fontsize_scalebar)
    x_kwargs = dict(verticalalignment='top',**common_font_inset)
    x_font,y_font = Scalebar.\
        font_kwargs_modified(x_kwargs=x_kwargs,
                             y_kwargs=dict(horizontalalignment='right',
                                           **common_font_inset))
    # set up the font, offset ('fudge') the text from the lines              
    fudge_x = dict(x=0,y=-0.05)
    fudge_y = dict(x=0,y=0.01)
    Scalebar.crossed_x_and_y_relative(0.42,0.57,ax=axins,
                                      x_kwargs=dict(width=2,unit="nm",
                                                    font_kwargs=x_font,
                                                    fudge_text_pct=fudge_x,
                                                    **kw_common),
                                      y_kwargs=dict(height=8,unit='kcal/\nmol',
                                                    font_kwargs=y_font,
                                                    fudge_text_pct=fudge_y,
                                                    **kw_common))
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    color_box = 'rebeccapurple'           
    PlotUtilities.color_frame('rebeccapurple',ax=axins) 
    Annotations.add_rectangle(ax_correction,xlim_box,ylim_box,
                              edgecolor=color_box)
    ax_correction.set_xlim(xlim_fec)
    ax_energy = plt.subplot(3,1,3)    
    plot_landscape(data,xlim_fec,kw_landscape=kw['kw_landscape'],
                   plot_derivative=False,label_deltaG=" ")
    ax_energy.set_xlim(xlim_fec)                         
    setup_pedagogy_ticks(ax_energy,scale_bar_x,x_heat_kw,y_heat_kw,
                         offset_y=offset_y_pedagogy)
    # add in the equation notation
    strings,colors = [],[]
    labels = kwargs_labels()
    # add in the appropriate symbols 
    strings = ["$\Delta G$ = ",labels[0]," + ",labels[1]," - ",labels[2]]
    colors_labels = [c['color'] for c in kwargs_correction()]
    colors = ["k"] + [item for list in [[c,"k"] for c in colors_labels]
                      for item in list]
    x,y = Scalebar.x_and_y_to_abs(x_rel=0.08,y_rel=0.85,ax=ax_energy)        
    Annotations.rainbow_text(x,y,strings=strings,colors=colors,
                             ax=ax_energy,size=legend_font_size)
    PlotUtilities.legend(handlelength=0.5,loc=(0.03,0.8))      
    PlotUtilities.no_x_label(ax_energy)           
    ax_to_label = [ax for i,ax in enumerate(fig.axes) if i != 2]
    PlotUtilities.label_tom(fig=fig,axis_func=lambda *a,**kw: ax_to_label,
                            loc=(-0.05,1.02))
    PlotUtilities.save_png_and_svg(fig,out_name)  
    

def print_info(helical_data):
    # print out the full DeltaG 
    full = helical_data[0].generate_landscape_obj()
    full_mean = full.mean_landscape_kcal_per_mol
    full_std = full.std_landscape_kcal_per_mol
    n_aa = 160
    print("DeltaG for the entire landscape is {:.3g} +/- {:.2g} kcal/mol".\
          format(full_mean[-1],full_std[-1]))
    print("DeltaG/aa for the entire landscape is {:.3g} +/- {:.2g} kcal/mol".\
          format(full_mean[-1]/n_aa,full_std[-1]/n_aa))
    # Get the ED helical part we care about
    ed = helical_data[1].generate_landscape_obj()
    ext_nm = ed._extension_grid_nm
    ed_mean_delta_per_nm = ed.mean_delta_landscape_kcal_per_mol_per_nm
    ed_std_delta_per_nm = ed.std_delta_landscape_kcal_per_mol_per_nm
    # get the nm we care about
    x_of_interest = 6.5
    idx_of_interest = np.argmin(abs(ext_nm - x_of_interest))
    ed_std_delta_per_nm_of_interest = ed_mean_delta_per_nm[idx_of_interest]
    ed_mean_delta_per_nm_of_interest = ed_std_delta_per_nm[idx_of_interest]
    mean_ed = ed_std_delta_per_nm_of_interest
    # use SEM / STD?
    n_bootstraps = 500
    aa_per_nm = 3
    error_ed = ed_mean_delta_per_nm_of_interest/np.sqrt(n_bootstraps)
    print("The top of the ED helix has energy {:.2f} +/- {:.2g} kcal/mol/aa".\
          format(mean_ed/aa_per_nm,error_ed/aa_per_nm))
    # get the maximim stiffness, in pN/nm
    stiffness_kcal_per_mol_nm_sq = \
        np.gradient(ed_mean_delta_per_nm)/np.gradient(ext_nm)
    # 1kT = 0.593 kcal/mol. 
    k_kT_per_nm_sq = stiffness_kcal_per_mol_nm_sq/0.593
    # 1kT = 4.1 pN * nm
    k_pN_per_nm_sq = k_kT_per_nm_sq * 4.1
    max_stiffness_pN_per_nm = max(k_pN_per_nm_sq)
    print("Max stiffness is {:.1f} kcal/mol".\
          format(max(stiffness_kcal_per_mol_nm_sq)))
    print("Max stiffness is {:.1f} pN/nm".format(max_stiffness_pN_per_nm))
    print("The max ED energy per nm is {:.2f} kcal/mol".\
          format(max(ed_mean_delta_per_nm)))
    cb = helical_data[2].generate_landscape_obj()
    print("The min ED energy per nm is {:.2f}/mol".\
          format(min(ed_mean_delta_per_nm)))


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    np.random.seed(42)
    flickering_dir = "../../LargerDataset/"
    # XXX use the flickering dir for stuff
    cache_dir = flickering_dir 
    force_recalculation = False
    GenUtilities.ensureDirExists(flickering_dir)
    bin_size_meters = 0.4e-9
    # write down the areas we want to look at 
    adhesion_min = 17e-9
    ed_max = 32e-9
    cd_max = 48e-9
    a_max = 70e-9
    slice_area = GenerateLandscapes.slice_area
    kw = dict(min_v_m=adhesion_min)
    areas = [\
        slice_area([adhesion_min,a_max],"Full (no adhesion)",**kw),
        slice_area([adhesion_min,ed_max],"Helix ED",**kw),
        slice_area([ed_max,cd_max],"Helix CB",**kw),
        slice_area([cd_max,a_max],"Helix A",**kw),
        ]    
    for a in areas:
        a.set_num_bins_by_bin_in_meters(bin_size_meters)         
    # read in the data 
    f = GenerateLandscapes.get_cacheable_data
    data_to_analyze = CheckpointUtilities.\
        getCheckpoint("./cached_landscapes.pkl",f,
                      force_recalculation,areas,flickering_dir,bin_size_meters)
    # split into the data we care about
    helical_data = []
    for a in areas:
        tmp = copy.deepcopy(data_to_analyze)
        l = copy.deepcopy(tmp.landscape[0][0])
        l.q -= min(l.q)
        min_v,max_v = a.ext_bounds_nm_rel*1e-9
        slice_idx = np.where( (l.q >= min_v) & (l.q <= max_v))[0]
        assert slice_idx.size > 0
        sanit = lambda x: x[slice_idx].copy()        
        flat_list = [tmp_landscape for list_v in tmp.landscape 
                     for tmp_landscape in list_v]
        for l in flat_list:
            l.q = sanit(l.q)
            # zero the energy
            l.energy = sanit(l.energy)
            min_e = min(l.energy)
            l.energy -=  min_e
            l.A_z = sanit(l.A_z)
            l.A_z_dot = sanit(l.A_z_dot)
            l.one_minus_A_z_ddot_over_k = sanit(l.one_minus_A_z_ddot_over_k)   
        helical_data.append(tmp)
    print_info(helical_data)
    make_pedagogical_plot(helical_data[0],landscape_kwargs()[0])
    # make the SI / second derivative plot
    make_second_deriv_plot(helical_data[0],landscape_kwargs()[0]) 
    # make the heatmaps/energy landscape plots
    make_detalied_plots(helical_data,areas)
    # make the 'gallery' plots.
    make_gallery_plot(areas,helical_data)


    
if __name__ == "__main__":
    run()
