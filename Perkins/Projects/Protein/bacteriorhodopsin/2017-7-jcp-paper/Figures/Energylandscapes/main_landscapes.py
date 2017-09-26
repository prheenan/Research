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
sys.path.append("../")
import jcp_fig_util

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

class cacheable_data:
    def __init__(self,landscape,heatmap_data,heatmap_data_z):
        self.landscape = landscape
        self.heatmap_data = heatmap_data
        self.heatmap_data_z = heatmap_data_z

class slice_area:
    def __init__(self,ext_bounds,plot_title):
        self.ext_bounds = ext_bounds
        self.plot_title = plot_title
        self.n_bins = None
    def set_num_bins_by_bin_in_meters(self,bin_size_meters):
        self.n_bins = \
            int(np.ceil(abs(np.diff(self.ext_bounds))/bin_size_meters))
    @property
    def save_name(self):
        return self.plot_title.replace(" ","_") + ".pkl"
        
        
class landscape_data:
    def __init__(self,landscape_objs,kT=4.1e-21):
        self.landscape_objs = landscape_objs
        self.kT = 4.1e-21
    @property        
    def _extension_grid_nm(self):
        # get the extensions...
        ext_nm = self._extensions_nm        
        # get the absolute min and max
        min_all = np.max([min(e) for e in ext_nm])
        max_all = np.min([max(e) for e in ext_nm])
        # get the grid to interpolate onto
        grid = np.linspace(min_all,max_all,num=len(ext_nm[0]),endpoint=True)   
        return grid   
    def from_Joules_to_kcal_per_mol(self):
        return IWT_Util.kT_to_kcal_per_mol() * (1/self.kT)
    def _raw_uninterpolared_landscapes_kcal_per_mol(self,l):
        return l.G_0 * self.from_Joules_to_kcal_per_mol()
    def _grid_property(self,property_func):
        grid = self._extension_grid_nm
        # get the raw data
        raw_y = [property_func(o) for o in self.landscape_objs]
        return grid_interpolate_arrays(self._extensions_nm,raw_y,grid)
    @property                
    def _landscapes_kcal_per_mol(self):
        """
        get the landscapes in kcal per mol 
        """
        func = self._raw_uninterpolared_landscapes_kcal_per_mol
        return self._grid_property(func)
    @property                        
    def _extensions_nm(self):
        return [l.q * 1e9 for l in self.landscape_objs]
    def amino_acids_per_nm(self):
        return 3     
    @property
    def _delta_landscapes_kcal_per_mol_per_AA(self):
        
        dy_kcal_mol = [np.gradient(y,**self.mean_std_opt()) 
                       for y in self._landscapes_kcal_per_mol]
        ext_nm = self._extensions_nm                   
        dx_nm = [np.gradient(x) for x in ext_nm]
        dx_aa = [d * self.amino_acids_per_nm() for d in dx_nm]
        # get the individual gradients
        gradients_per_aa = [(dy_i/(dx_i))  
                            for dy_i,dx_i in zip(dy_kcal_mol,dx_aa)]
        grid_x = self._extension_grid_nm                            
        to_ret =  grid_interpolate_arrays(ext_nm,gradients_per_aa,grid_x)
        return to_ret 
    def mean_std_opt(self):
        return dict(axis=0)
    @property        
    def mean_landscape_kcal_per_mol(self):
        # get the landscapes, XXX need to interpolate back onto uniform grid. 
        return np.mean(self._landscapes_kcal_per_mol,**self.mean_std_opt())
    @property                
    def std_landscape_kcal_per_mol(self):
        return np.std(self._landscapes_kcal_per_mol,**self.mean_std_opt())
    @property        
    def mean_delta_landscape_kcal_per_mol_per_AA(self):
        return np.mean(self._delta_landscapes_kcal_per_mol_per_AA,
                       **self.mean_std_opt())
    @property                               
    def std_delta_landscape_kcal_per_mol_per_AA(self):
        return np.std(self._delta_landscapes_kcal_per_mol_per_AA,
                       **self.mean_std_opt())     

def grid_interpolate_arrays(x_arr,y_arr,x_grid):
    to_ret = []
    for x, y in zip(x_arr,y_arr):
        # get the interpolated data
        tmp_grid = scipy.interpolate.griddata(points=x,
                                              values=y,
                                              xi=x_grid)
        to_ret.append(tmp_grid)                                     
    return to_ret
                       
def heatmap(x,y,bins):                       
    # concatenate everything
    cat_x = np.concatenate(x)
    cat_y = np.concatenate(y)
    # SEE: histogram2d documentation
    histogram, x_edges,y_edges = \
        np.histogram2d(cat_x,cat_y,bins=bins)
    # each row should list y; transpose so this is the case 
    histogram = histogram.T            
    return histogram, x_edges,y_edges 
    
def get_heatmap_data(time_sep_force_arr,bins=(100,100)):        
    sep_nm = [t.Separation*1e9 for t in time_sep_force_arr]
    z_nm = [t.ZSnsr*1e9 for t in time_sep_force_arr]
    force_pN = [t.Force*1e12 for t in time_sep_force_arr]
    heatmap_force_extension = heatmap(sep_nm,force_pN,bins=bins)
    heatmap_force_z = heatmap(z_nm,force_pN,bins=bins)
    return heatmap_force_extension,heatmap_force_z

    
def _get_landscapes(iwt_obj_subsets,n_bins):    
    for objs in iwt_obj_subsets:
        yield InverseWeierstrass.free_energy_inverse_weierstrass(objs)
            
def get_area_bounds(objs,area):
    z_0,z_1 = area.ext_bounds
    average_v = np.mean([r.Velocity for r in objs])    
    dt_step = objs[0].Time[1] - objs[0].Time[0]
    z_0_arr = [ (np.where(o.Separation <= z_0)[0][-1]) for o in objs]
    average_idx_delta = int(np.round((z_1-z_0)/(average_v *dt_step)))
    # determine how large the delta can actually, so all the objects
    # lie on the grid
    sizes = [o.Force.size for o in objs]
    actual_delta = min([min(average_idx_delta,s-z_tmp-1)
                        for z_tmp,s in zip(z_0_arr,sizes)])
    z_f_arr = [z + actual_delta for z in z_0_arr]
    to_ret = [ [z0,zf] for z0,zf in zip(z_0_arr,z_f_arr)]
    for i,bounds in enumerate(to_ret):
        np.testing.assert_allclose(np.diff(bounds),np.diff(to_ret[0]))
        n = objs[i].Force.size
        assert(bounds[-1] < n) , \
            "{:d}/{:d}".format(bounds[-1],n)
    # POST: all the bounds match 
    return to_ret
    

    
def filter_landscapes(landscapes,n_bins):
    # zero everything; don't care about absolute offsets.
    for l in landscapes:
        l.q -= min(l.q)
    min_v = 0
    max_v = max([max(l.q) for l in landscapes])
    bins = np.linspace(min_v,max_v,n_bins)
    to_ret = [WeierstrassUtil._filter_single_landscape(l,bins) 
              for l in landscapes]
    return to_ret

def get_cacheable_data(areas,flickering_dir,heat_bins=(100,100),
                       offset_N=7.1e-12):
    raw_data = IoUtilHao.read_and_cache_data_hao(None,force=False,
                                                 cache_directory=flickering_dir,
                                                 limit=10,
                                                 renormalize=False)
    # only look at data with ~300nm/s
    v_exp = 300e-9
    raw_data = [r for r in raw_data 
                if np.allclose(r.Velocity,v_exp,atol=0,rtol=0.01)]
    raw_area_slices = [[] for _ in areas]
    area_bounds = [get_area_bounds(raw_data,area) for area in areas]
    # fix all the spring constants. XXX need to account for this...
    mean_spring_constant = np.mean([r.LowResData.meta.SpringConstant for r in raw_data])
    for r in raw_data:
        r.LowResData.meta.SpringConstant = mean_spring_constant
    for i,r in enumerate(raw_data):
        r.Force -= offset_N
        for j,area in enumerate(areas):
           idx_0,idx_f = area_bounds[j][i]
           s = slice(idx_0,idx_f,1)
           this_area = FEC_Util.MakeTimeSepForceFromSlice(r,s)
           raw_area_slices[j].append(this_area)
        # r is no longer needed; stop referencing it to make space
        raw_data[i] = None
    to_ret = []
    skip = 0
    N_boostraps = 3
    for area,slice_tmp in zip(areas,raw_area_slices):
        # for each area, use the same starting seed 
        # (that the data are consistent)
        np.random.seed(42)
        # get the heatmap histograms
        heatmap_data = get_heatmap_data(slice_tmp)  
        # get the landscapes (we use N, to get an error)
        iwt_objs = WeierstrassUtil.convert_list_to_iwt(slice_tmp)
        # delete the original list to free its memnory
        slice_tmp[:] = []
        n_objs = len(iwt_objs)
        ids = np.arange(n_objs,dtype=np.int64)
        # randomly choose the ids with replacement for bootstrapping\
        choose_ids = lambda : np.random.choice(ids,size=n_objs,replace=True)
        # skip the first N (if we already chose those, assuming a consistent 
        # seed)
        skipped = [choose_ids() for i in range(skip)]
        # get the actual ids we want 
        id_choices = [choose_ids() for i in range(N_boostraps)]
        iwt_obj_subsets = [ [iwt_objs[i] for i in a] for a in id_choices]
        functor = lambda : _get_landscapes(iwt_obj_subsets,area.n_bins)
        name_func = lambda  i,d: \
            area.save_name + "_bootstrap_{:d}".format(i+skip) 
        cache_dir = "./{:s}/".format(area.save_name)
        GenUtilities.ensureDirExists(cache_dir)
        iwt_tmp = CheckpointUtilities.multi_load(cache_dir,load_func=functor,
                                                 force=False,
                                                 name_func=name_func)
        # get the filtered version of the landscape, using the specified bin
        # size 
        filtered_iwt = filter_landscapes(iwt_tmp,area.n_bins)
        # make the object we want for this 'area' slice
        to_ret.append(cacheable_data(filtered_iwt,*heatmap_data))
    return to_ret
    
def make_heatmap(histogram, x_edges,y_edges,kw_heatmap):
    # XXX ? digitize all the ids so we know what bin they fall into...
    X,Y = np.meshgrid(x_edges,y_edges)
    X -= np.min(X)
    plt.gca().pcolormesh(X,Y,histogram,**kw_heatmap)
        
def plot_landscape(data,xlim,kw_landscape=dict(),plot_derivative=True,
                   label_deltaG = PlotUtilities.variable_string("\Delta G")):
    landscape_kcal_per_mol = data.mean_landscape_kcal_per_mol
    landscape_kcal_per_mol -= min(landscape_kcal_per_mol)
    std_landscape_kcal_per_mol = data.std_landscape_kcal_per_mol
    extension_nm = data._extension_grid_nm
    extension_nm -= min(extension_nm)
    extension_aa = data.amino_acids_per_nm() * extension_nm
    grad = lambda x: np.gradient(x)/(np.gradient(extension_aa))
    delta_landscape_kcal_per_mol_per_amino_acid = grad(landscape_kcal_per_mol)
    landscape_upper = landscape_kcal_per_mol+std_landscape_kcal_per_mol
    landscape_lower =landscape_kcal_per_mol-std_landscape_kcal_per_mol
    delta_landscape_kcal_per_mol_per_amino_acid = \
        data.mean_delta_landscape_kcal_per_mol_per_AA
    std_delta_landscape_kcal_per_mol_per_AA = \
        data.std_delta_landscape_kcal_per_mol_per_AA
    upper_delta_landscape = delta_landscape_kcal_per_mol_per_amino_acid+\
                            std_delta_landscape_kcal_per_mol_per_AA
    lower_delta_landscape = delta_landscape_kcal_per_mol_per_amino_acid-\
                            std_delta_landscape_kcal_per_mol_per_AA
    # make a second axis for the number of ammino acids 
    limits_delta = [min(delta_landscape_kcal_per_mol_per_amino_acid),
                    max(delta_landscape_kcal_per_mol_per_amino_acid)]    
    limits_energy = [min(landscape_kcal_per_mol),max(landscape_kcal_per_mol)]
    units_energy_delta = label_deltaG + \
                         r" per AA (kcal/(mol $\cdot$ AA))"   
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
                                        label=units_energy,color=landscape_color,
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
                          xlim=None): 
    """
    Creates a plot of
    """
    heatmap_data = data_to_plot.heatmap_data
    data_landscape = landscape_data(data_to_plot.landscape)
        
    # # ploy the heat map 
    ax_heat = plt.subplot(2,1,1)
    heatmap_plot(heatmap_data,data_landscape.amino_acids_per_nm(),
                 kw_heatmap=kw_heatmap)
    if (xlim is None):                 
        xlim = plt.xlim()
    # # plot the energy landscape...
    ax_energy = plt.subplot(2,1,2)    
    plot_landscape(data_landscape,xlim,kw_landscape=kw_landscape)
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
    for i,d in enumerate(data_to_analyze):
        fig = PlotUtilities.figure((3.25,5))     
        ax = create_landscape_plot(d,xlim=[None,None],**(kwargs[i]))
        out_name = "landscape{:d}_{:s}".format(i,areas[i].plot_title)
        axis_func = lambda x: [x[0],x[2]]
        PlotUtilities.label_tom(fig,axis_func=axis_func)
        PlotUtilities.save_png_and_svg(fig,out_name.replace(" ","_"))    
    
def normalize_axes(ax_list,manual_min=None):
    """
    ensures all the y axes in ax_list have the same limits
    """
    max_ylim = np.max([ax.get_ylim() for ax in ax_list])
    min_ylim = np.min([ax.get_ylim() for ax in ax_list])
    range = max_ylim-min_ylim
    fudge = range * 0.02
    min_v = min_ylim-fudge
    if (manual_min is not None):
        min_v = manual_min    
    for ax in ax_list:
        ax.set_ylim([min_v,max_ylim+fudge])    
    
def helical_gallery_plot(helical_areas,helical_data,helical_kwargs):
    axs,first_axs,second_axs = [],[],[]
    offset_y = 0.05
    kw_scalebars = [dict(offset_x=0.5,offset_y=offset_y),
                    dict(offset_x=0.35,offset_y=offset_y),
                    dict(offset_x=0.5,offset_y=offset_y)]
    xlims = [ [None,None],[None,None],[None,59.1]    ]   
    arrow_x = [0.55,0.5,0.60]
    arrow_y = [0.60,0.6,0.5]
    for i,(a,data) in enumerate(zip(helical_areas,helical_data)):
        kw_tmp = helical_kwargs[i]
        data_landscape = landscape_data(data.landscape)
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
        ax_1.annotate("",xytext=(arrow_x[i],arrow_y[i]),textcoords='axes fraction',
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
        PlotUtilities.title(a.plot_title,color=color)
    # make all the axes consistent            
    normalize_axes(first_axs)
    normalize_axes(second_axs)
    # determine where the min need to be for the zeros to be OK
    min_y,max_y = first_axs[0].get_ylim()
    normalized_zero = -min_y/(max_y-min_y)
    min_y_second,max_y_second = second_axs[0].get_ylim()
    # want f = 0-min/(max-min) --> min = f*max/(f-1)
    needed_min = normalized_zero * max_y_second / (normalized_zero-1)
    normalize_axes(second_axs,manual_min=needed_min)    
    # after normalization, add in the scale bars 
    for i,(ax_1,ax_2) in enumerate(zip(first_axs,second_axs)):
        Scalebar.x_scale_bar_and_ticks_relative(unit="nm",width=5,ax=ax_2,
                                                **kw_scalebars[i])
        PlotUtilities.no_x_label(ax_2)     
            
    
    
def make_gallery_plot(areas,data_to_analyze,out_name="./gallery"):
    # skip the first one (the entire landscape )
    helical_areas = areas[1:]
    helical_data = data_to_analyze[1:]
    helical_kwargs = landscape_kwargs()[1:]
    fig = PlotUtilities.figure((7,2.5))
    helical_gallery_plot(helical_areas,helical_data,helical_kwargs)    
    PlotUtilities.save_png_and_svg(fig,out_name)                    
    
def setup_pedagogy_ticks(ax,scale_bar_x,x_heat_kw,y_heat_kw,offset_y=0.9):
    font_kwargs= copy.deepcopy(Scalebar.def_font_kwargs_y)
    # fix the keywords relative to the heatmap
    x_heat_kw['font_kwargs']['color'] = 'k'
    x_heat_kw['line_kwargs']['color'] = 'k'
    y_heat_kw['font_kwargs']['color'] = 'k'
    y_heat_kw['line_kwargs']['color'] = 'k'  
    y_heat_kw['unit'] = 'kcal/mol '
    y_heat_kw['height'] = 40
    Scalebar.crossed_x_and_y_relative(scale_bar_x,offset_y,ax=ax,
                                      x_kwargs=x_heat_kw,
                                      y_kwargs=y_heat_kw)
    PlotUtilities.no_y_label(ax)                                              

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
    ext_nm = data._extension_grid_nm.copy()
    ext_nm -= min(ext_nm)
    convert = data.from_Joules_to_kcal_per_mol()
    energies = [data._grid_property(lambda x: x.A_z * convert),
                data._grid_property(lambda x: -1 * x.first_deriv_term * convert),
                data._grid_property(lambda x: x.second_deriv_term* convert)]
    labels = kwargs_labels()
    kwargs = kwargs_correction()
    landscape_kcal_per_mol = data.mean_landscape_kcal_per_mol             
    for i,e in enumerate(energies):
        energy_rel = np.mean(e,axis=0)
        plt.plot(ext_nm,energy_rel-min(energy_rel),label=labels[i],**kwargs[i])  
    
def make_pedagogical_plot(data_to_plot,kw,out_name="./iwt_diagram"):
    heatmap_data = data_to_plot.heatmap_data
    data = landscape_data(data_to_plot.landscape)
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
    common_kw = dict(color='w')
    x_font,y_font = Scalebar.\
        font_kwargs_modified(x_kwargs=common_kw,
                             y_kwargs=common_kw)
    heat_kw_common = dict(line_kwargs=dict(color='w',linewidth=1.5))
    x_heat_kw = dict(width=15,unit="nm",font_kwargs=x_font,**heat_kw_common)
    y_heat_kw = dict(height=30,unit='pN ',font_kwargs=y_font,**heat_kw_common)
    # add a scale bar for the heatmap...
    scale_bar_x = 0.83
    Scalebar.crossed_x_and_y_relative(scale_bar_x,0.55,ax=ax_heat,
                                      x_kwargs=x_heat_kw,
                                      y_kwargs=y_heat_kw)
    jcp_fig_util.add_helical_boxes(ax=ax_heat,ymax_box=0.9,alpha=1.0,
                                   font_color='w')
    # # plot the energy landscape...
    ax_correction = plt.subplot(3,1,2)    
    plot_with_corrections(data)
    PlotUtilities.no_x_label(ax_correction)
    PlotUtilities.lazyLabel("","Energy (kcal/mol)","")
    ax_correction.set_xlim(xlim_fec)            
    offset_y_pedagogy = 0.40
    setup_pedagogy_ticks(ax_correction,scale_bar_x,x_heat_kw,y_heat_kw,
                         offset_y=offset_y_pedagogy)
    legend_font_size = 9                         
    legend = PlotUtilities.legend(handlelength=1.5,loc=(0.15,0.07),ncol=3,
                                  fontsize=legend_font_size,handletextpad=0.4)
    for i,text in enumerate(legend.get_texts()):
        plt.setp(text, color = kwargs_correction()[i]['color'])    
    # make the inset plot 
    axins = zoomed_inset_axes(ax_correction, zoom=4, loc=2,
                              borderpad=0.8) 
    plot_with_corrections(data)
    xlim_box = [1,6]
    ylim_box = [-3,15]
    plt.xlim(xlim_box)
    plt.ylim(ylim_box)
    PlotUtilities.no_x_anything(axins)
    PlotUtilities.no_y_anything(axins)
    # add in scale bars
    kw_common = dict(line_kwargs=dict(linewidth=0.75,color='k'))
    common_font_inset = dict(fontsize=6)
    x_kwargs = dict(verticalalignment='top',**common_font_inset)
    x_font,y_font = Scalebar.\
        font_kwargs_modified(x_kwargs=x_kwargs,
                             y_kwargs=dict(horizontalalignment='left',
                                           **common_font_inset))
    # set up the font, offset ('fudge') the text from the lines                              
    fudge_x = dict(x=0,y=-0.5)
    fudge_y = dict(x=0,y=0.1)
    Scalebar.crossed_x_and_y_relative(0.82,0.48,ax=axins,
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
    Annotations.add_rectangle(ax_correction,xlim_box,ylim_box,edgecolor=color_box)
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
    PlotUtilities.save_png_and_svg(fig,out_name)  
    
    
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
    bin_size_meters = 0.2e-9
    # write down the areas we want to look at 
    adhesion_min = 17e-9
    ed_max = 32e-9
    cd_max = 48e-9
    a_max = 65e-9
    areas = [\
        slice_area([adhesion_min,a_max],"Full (no adhesion)"),
        slice_area([adhesion_min,ed_max],"Helix ED"),
        slice_area([ed_max,cd_max],"Helix CB"),
        slice_area([cd_max,a_max],"Helix A"),
        ]    
    for a in areas:
        a.set_num_bins_by_bin_in_meters(bin_size_meters)         
    # read in the data 
    data_to_analyze = CheckpointUtilities.\
        getCheckpoint("./cached_landscapes.pkl",get_cacheable_data,
                      force_recalculation,areas,flickering_dir,bin_size_meters)
    make_pedagogical_plot(data_to_analyze[0],landscape_kwargs()[0])
    # make the heatmaps/energy landscape plots
    make_detalied_plots(data_to_analyze,areas)
    # make the 'gallery' plots.
    make_gallery_plot(areas,data_to_analyze)
    
if __name__ == "__main__":
    run()
