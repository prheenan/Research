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
import scipy

class cacheable_data:
    def __init__(self,landscape,heatmap_data):
        self.landscape = landscape
        self.heatmap_data = heatmap_data

class slice_area:
    def __init__(self,ext_bounds,plot_title,n_bins):
        self.ext_bounds = ext_bounds
        self.plot_title = plot_title
        self.n_bins = n_bins
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
        assert sum([len(e) == len(ext_nm[0]) for e in ext_nm]) == len(ext_nm) ,\
            "Not all landscapes binned the same way..."            
        # get the absolute min and max
        min_all = np.max([min(e) for e in ext_nm])
        max_all = np.min([max(e) for e in ext_nm])
        # get the grid to interpolate onto
        grid = np.linspace(min_all,max_all,num=len(ext_nm[0]),endpoint=True)   
        return grid   
    @property    
    def _raw_uninterpolared_landscapes_kcal_per_mol(self):
        return [IWT_Util.kT_to_kcal_per_mol() * (l.EnergyLandscape/self.kT)
                for l in self.landscape_objs]
    @property                
    def _landscapes_kcal_per_mol(self):
        grid = self._extension_grid_nm
        # get the raw data
        raw_y = self._raw_uninterpolared_landscapes_kcal_per_mol
        return grid_interpolate_arrays(self._extensions_nm,raw_y,grid)
    @property                        
    def _extensions_nm(self):
        return [l.Extensions * 1e9 for l in self.landscape_objs]
    def amino_acids_per_nm(self):
        return 3     
    @property
    def _delta_landscapes_kcal_per_mol_per_AA(self):
        dy = np.gradient(self._raw_uninterpolared_landscapes_kcal_per_mol,
                         **self.mean_std_opt())
        ext_nm = self._extensions_nm                   
        dx = np.gradient(ext_nm,**self.mean_std_opt())
        # get the individual gradients
        gradients_per_aa = [(dy_i/(dx_i))  * self.amino_acids_per_nm()
                            for dy_i,dx_i in zip(dy,dx)]
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
                       
def get_heatmap_data(time_sep_force_arr,bins=(100,100)):        
    sep_nm = [t.Separation*1e9 for t in time_sep_force_arr]
    force_pN = [t.Force*1e12 for t in time_sep_force_arr]
    id_array = [ i for i,_ in enumerate(time_sep_force_arr)]
    # concatenate everything
    cat_sep_nm = np.concatenate(sep_nm)
    cat_force_pN = np.concatenate(force_pN)
    # SEE: histogram2d documentation
    histogram, x_edges,y_edges = \
        np.histogram2d(cat_sep_nm,cat_force_pN,bins=bins)
    # each row should list y; transpose so this is the case 
    histogram = histogram.T            
    return histogram, x_edges,y_edges 
    
def get_cacheable_data(areas,flickering_dir,heat_bins=(100,100),
                       offset_N=7.1e-12):
    force_read_data = False    
    raw_data = IoUtilHao.read_and_cache_data_hao(None,force=force_read_data,
                                                 cache_directory=flickering_dir,
                                                 limit=None,
                                                 renormalize=False)
    raw_area_slices = [[] for _ in areas]
    for i,r in enumerate(raw_data):
        r.Force -= offset_N
        for j,area in enumerate(areas):
           this_area = FEC_Util.slice_by_separation(r,*area.ext_bounds)
           raw_area_slices[j].append(this_area)
        # r is no longer needed; stop referencing it to make space
        raw_data[i] = None
    to_ret = []
    N = 3
    for area,slice_tmp in zip(areas,raw_area_slices):
        # get the heatmap histograms
        heatmap_data = get_heatmap_data(slice_tmp)  
        # get the landscapes (we use N, to get an error)
        iwt_objs = IWT_Util.convert_list_to_iwt(slice_tmp)
        # delete the original list to free its memnory
        slice_tmp[:] = []
        n_objs = len(iwt_objs)
        num_n = int(np.round(n_objs/N))
        assert (n_objs % N == 0), \
            "Need a multiple of {:d}, got {:d}".format(N,n_objs)
        ids = [i for i in range(len(iwt_objs))]
        # shuffle ids in-place
        np.random.shuffle(ids)
        # get the groups; 
        # (1) ids was 0,1,2,...
        # (2) it is shuffled to (e.g.) 8,1,3,0
        # (3) we want to group (if N=2), like [ [8,1],[3,0], ...]
        iwt_obj_subsets = [ [iwt_objs[j] 
                             for j in ids[N*i:N*i+N]]
                            for i in range(0,num_n)]
        iwt_helix_data_tmp = \
            [InverseWeierstrass.FreeEnergyAtZeroForce(objs,area.n_bins)
             for objs in iwt_obj_subsets]
        # make the object we want for this 'area' slice
        to_ret.append(cacheable_data(iwt_helix_data_tmp,heatmap_data))
    return to_ret
    
def make_heatmap(histogram, x_edges,y_edges,kw_heatmap):
    # XXX ? digitize all the ids so we know what bin they fall into...
    X,Y = np.meshgrid(x_edges,y_edges)
    plt.gca().pcolormesh(X,Y,histogram,**kw_heatmap)
    
def plot_landscape(data,xlim,kw_landscape=dict()):
    landscape_kcal_per_mol = data.mean_landscape_kcal_per_mol
    std_landscape_kcal_per_mol = data.std_landscape_kcal_per_mol
    extension_nm = data._extension_grid_nm
    extension_aa = data.amino_acids_per_nm() * extension_nm
    grad = lambda x: np.gradient(x)/(np.gradient(extension_aa))
    delta_landscape_kcal_per_mol_per_amino_acid = grad(landscape_kcal_per_mol)
    landscape_upper = landscape_kcal_per_mol+std_landscape_kcal_per_mol
    landscape_lower =landscape_kcal_per_mol-std_landscape_kcal_per_mol 
    upper_delta_landscape = grad(landscape_upper)
    lower_delta_landscape = grad(landscape_lower)
    ax_energy = plt.gca()
    # plot the landscape and its standard deviation
    plt.plot(extension_nm,landscape_kcal_per_mol,
             linestyle='--',**kw_landscape)
    plt.fill_between(x=extension_nm,
                     y1=landscape_lower,
                     y2=landscape_upper,
                     alpha=0.3,**kw_landscape)
    # make a second axis for the number of ammino acids 
    units_energy = r"($\frac{\mathrm{kcal}}{\mathrm{mol}}$)"
    units_energy_delta = r"($\frac{\mathrm{kcal}}{\mathrm{mol} \cdot AA}$)"
    PlotUtilities.lazyLabel("Extension (nm)","Free energy " + units_energy,"")    
    limits_delta = [min(delta_landscape_kcal_per_mol_per_amino_acid),
                    max(delta_landscape_kcal_per_mol_per_amino_acid)]
    label = "Free energy difference " + units_energy_delta
    difference_color = 'rebeccapurple'
    ax_2 = PlotUtilities.secondAxis(ax_energy,
                                    label=label,color=difference_color,
                                    limits=limits_delta,secondY =True)
    # plot the energy delta and its bounds, based on the bounds on the landscape    
    ax_2.plot(extension_nm,delta_landscape_kcal_per_mol_per_amino_acid,
              color=difference_color,linestyle='-',linewidth=0.5)    
              
def heatmap_plot(heatmap_data,amino_acids_per_nm,kw_heatmap=dict()):
    ax_heat = plt.gca()
    make_heatmap(*heatmap_data,kw_heatmap=kw_heatmap)
    PlotUtilities.lazyLabel("","Force (pN)","")
    # make a second x axis for the number of ammino acids 
    xlim_fec = plt.xlim()
    limits = np.array(xlim_fec) * amino_acids_per_nm
    PlotUtilities.secondAxis(ax_heat,"Extension (AA #)",limits,secondY =False)
    plt.xlim(xlim_fec)
    PlotUtilities.no_x_label(ax_heat)    
    
def create_landscape_plot(data_to_plot,kw_heatmap=dict(),kw_landscape=dict()): 
    heatmap_data = data_to_plot.heatmap_data
    data_landscape = landscape_data(data_to_plot.landscape)
    # # ploy the heat map 
    ax_heat = plt.subplot(2,1,1)
    heatmap_plot(heatmap_data,data_landscape.amino_acids_per_nm(),
                 kw_heatmap=kw_heatmap)
    xlim_fec = plt.xlim()
    # # plot the energy landscape...
    ax_energy = plt.subplot(2,1,2)    
    plot_landscape(data_landscape,xlim_fec,kw_landscape=kw_landscape)
    
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
    n_bins = 250
    n_bins_helix_a = 100
    n_bins_helix_e = 100
    # write down the areas we want to look at 
    areas = [\
        slice_area([18e-9,75e-9],"Full (no adhesion)",n_bins),
        slice_area([18e-9,27e-9],"Helix E",n_bins_helix_e),
        slice_area([50e-9,75e-9],"Helix A",n_bins_helix_a),
        ]    
    # read in the data 
    data_to_analyze = CheckpointUtilities.\
        getCheckpoint("./cached_landscapes.pkl",get_cacheable_data,
                      force_recalculation,areas,flickering_dir)
    kwargs = [ dict(kw_heatmap=dict(cmap=plt.cm.Greys_r),
                    kw_landscape=dict(color='k')),
               dict(kw_heatmap=dict(cmap=plt.cm.Blues_r),
                    kw_landscape=dict(color='royalblue')),
               dict(kw_heatmap=dict(cmap=plt.cm.Greens_r),
                    kw_landscape=dict(color='g'))]
    for i,d in enumerate(data_to_analyze):
        fig = PlotUtilities.figure((3.25,7))     
        create_landscape_plot(d,**(kwargs[i]))
        out_name = "landscape{:d}_{:s}".format(i,areas[i].plot_title)
        axis_func = lambda x: [x[0],x[2]]
        PlotUtilities.label_tom(fig,axis_func=axis_func)
        PlotUtilities.save_png_and_svg(fig,out_name.replace(" ","_"))
        

    
if __name__ == "__main__":
    run()
