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
    def nanometers_to_amino_acids(self):
        return 0.3        
    @property
    def _delta_landscapes_kcal_per_mol_per_AA(self):
        dy = np.gradient(self._raw_uninterpolared_landscapes_kcal_per_mol,
                         **self.mean_std_opt())
        ext_nm = self._extensions_nm                   
        dx = np.gradient(ext_nm,**self.mean_std_opt())
        # get the individual gradients
        gradients_per_aa = [dy_i/(self.nanometers_to_amino_acids()*dx_i) 
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
    
def get_cacheable_data(areas,flickering_dir,heat_bins=(100,100)):
    force_read_data = False    
    raw_data = IoUtilHao.read_and_cache_data_hao(None,force=force_read_data,
                                                 cache_directory=flickering_dir,
                                                 limit=None,renormalize=True)
    raw_area_slices = []
    for area in areas:
        this_area = [FEC_Util.slice_by_separation(r,*area.ext_bounds) 
                     for r in raw_data]
        raw_area_slices.append(this_area)
    to_ret = []
    N = 3
    for area,slice_tmp in zip(areas,raw_area_slices):
        # get the heatmap histograms
        heatmap_data = get_heatmap_data(slice_tmp,)  
        # get the landscapes (we use N, to get an error)
        iwt_objs = IWT_Util.convert_to_iwt(slice_tmp)
        num_n = int(np.round(len(iwt_objs)/N))
        assert (num_n % N == 0), \
            "Need a multiple of {:d}".format(N)
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
    
def make_heatmap(histogram, x_edges,y_edges):
    # XXX ? digitize all the ids so we know what bin they fall into...
    X,Y = np.meshgrid(x_edges,y_edges)
    plt.gca().pcolormesh(X,Y,histogram,cmap=plt.cm.afmhot)
    
def plot_landscape(data,xlim):
    landscape_kcal_per_mol = data.mean_landscape_kcal_per_mol
    extension_nm = data._extension_grid_nm
    delta_landscape_kcal_per_mol_per_amino_acid = \
        data.mean_delta_landscape_kcal_per_mol_per_AA
    std_landscape_kcal_per_mol = data.std_landscape_kcal_per_mol
    std_delta_landscape = data.std_delta_landscape_kcal_per_mol_per_AA
    ax_energy = plt.gca()
    plt.plot(extension_nm,landscape_kcal_per_mol,color='k')
    plt.fill_between(x=extension_nm,
                     y1=landscape_kcal_per_mol-std_landscape_kcal_per_mol,
                     y2=landscape_kcal_per_mol+std_landscape_kcal_per_mol,
                     color='k',alpha=0.3)
    # make a second axis for the number of ammino acids 
    units_energy = r"($\frac{\mathrm{kcal}}{\mathrm{mol}}$)"
    units_energy_delta = r"($\frac{\mathrm{kcal}}{\mathrm{mol} \cdot AA}$)"
    PlotUtilities.lazyLabel("Extension (nm)","Free energy " + units_energy,"")    
    limits_delta = [min(delta_landscape_kcal_per_mol_per_amino_acid),
                    max(delta_landscape_kcal_per_mol_per_amino_acid)]
    label = "Free energy difference " + units_energy_delta
    ax_2 = PlotUtilities.secondAxis(ax_energy,
                                    label=label,color='r',
                                    limits=limits_delta,secondY =True)
    ax_2.plot(extension_nm,delta_landscape_kcal_per_mol_per_amino_acid,
              color='r',linestyle='-',linewidth=0.5)                               

def heatmap_plot(heatmap_data,nanometers_to_amino_acids):
    ax_heat = plt.gca()
    make_heatmap(*heatmap_data)
    PlotUtilities.lazyLabel("","Force (pN)","")
    # make a second x axis for the number of ammino acids 
    xlim_fec = plt.xlim()
    limits = np.array(xlim_fec) * nanometers_to_amino_acids
    PlotUtilities.secondAxis(ax_heat,"Extension (AA #)",limits,secondY =False)
    plt.xlim(xlim_fec)
    PlotUtilities.no_x_label(ax_heat)    
    
def create_landscape_plot(data_to_plot): 
    heatmap_data = data_to_plot.heatmap_data
    data_landscape = landscape_data(data_to_plot.landscape)
    """
    for x,d in zip(data_landscape._extensions_nm,
                   data_landscape._raw_uninterpolared_landscapes_kcal_per_mol()):
        plt.plot(x,d)
    plt.plot(data_landscape._extension_grid_nm,
             data_landscape.mean_landscape_kcal_per_mol,
             'r--')
    one_std_lower = data_landscape.mean_landscape_kcal_per_mol - \
                    data_landscape.std_landscape_kcal_per_mol
    one_std_higher = data_landscape.mean_landscape_kcal_per_mol + \
                    data_landscape.std_landscape_kcal_per_mol
    plt.fill_between(x=data_landscape._extension_grid_nm,
                     y1=one_std_lower,y2=one_std_higher,color='k',alpha=0.3)           
    extension_nm = data._extension_grid_nm
    delta_landscape_kcal_per_mol_per_amino_acid = \
        data.mean_delta_landscape_kcal_per_mol_per_AA
    std_landscape_kcal_per_mol = data.std_landscape_kcal_per_mol
    std_delta_landscape = data.std_delta_landscape_kcal_per_mol_per_AA        
    plt.show()
    """
    
    # # ploy the heat map 
    ax_heat = plt.subplot(2,1,1)
    heatmap_plot(heatmap_data,data_landscape.nanometers_to_amino_acids())
    xlim_fec = plt.xlim()
    # # plot the energy landscape...
    ax_energy = plt.subplot(2,1,2)    
    plot_landscape(data_landscape,xlim_fec)
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    np.random.seed(42)
    flickering_dir = "../Data/"
    # XXX use the flickering dir for stuff
    cache_dir = flickering_dir 
    force_recalculation = False
    GenUtilities.ensureDirExists(flickering_dir)
    n_bins = 150
    n_bins_helix_a = 50
    n_binx_helix_e = 75
    # write down the areas we want to look at 
    areas = [\
        slice_area([18e-9,75e-9],"Full (no adhesion)",n_bins),
        slice_area([20e-9,27e-9],"Helix A",n_bins_helix_a),
        slice_area([50e-9,75e-9],"Helix E",n_binx_helix_e),
        ]    
    # read in the data 
    data_to_analyze = CheckpointUtilities.\
        getCheckpoint("./cached_landscapes.pkl",get_cacheable_data,
                      force_recalculation,areas,flickering_dir)
    for i,d in enumerate(data_to_analyze):
        fig = PlotUtilities.figure((3.25,7))        
        create_landscape_plot(d)
        out_name = "landscape{:d}_{:s}.png".format(i,areas[i].plot_title)
        PlotUtilities.savefig(fig,out_name)

    
if __name__ == "__main__":
    run()
