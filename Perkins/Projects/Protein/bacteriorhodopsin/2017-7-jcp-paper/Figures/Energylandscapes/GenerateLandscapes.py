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

import scipy,copy
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao


from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass,WeierstrassUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot    
from GeneralUtil.python import CheckpointUtilities,PlotUtilities,GenUtilities

class cacheable_data(object):
    def __init__(self,landscape,heatmap_data,heatmap_data_z):
        self.landscape = landscape
        self.heatmap_data = heatmap_data
        self.heatmap_data_z = heatmap_data_z
    def __deepcopy__(self, memo={}):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k,v in self.__dict__.items():
            setattr(result,k,copy.deepcopy(v,memo))
        return result

class slice_area(object):
    def __init__(self,ext_bounds,plot_title,min_v_m):
        self.ext_bounds = np.array(ext_bounds)
        self.plot_title = plot_title
        self.n_bins = None
        self.min_v_m = min_v_m
    def set_num_bins_by_bin_in_meters(self,bin_size_meters):
        self.n_bins = \
            int(np.ceil(abs(np.diff(self.ext_bounds))/bin_size_meters))
    @property
    def ext_bounds_nm(self):
        return self.ext_bounds * 1e9
    @property
    def ext_bounds_nm_rel(self):
        return self.ext_bounds_nm - (self.min_v_m*1e9)
    @property
    def save_name(self):
        return self.plot_title.replace(" ","_") + ".pkl"
        
        
class landscape_data(object):
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
    max_v = min([max(l.q) for l in landscapes])
    bins = np.linspace(min_v,max_v,n_bins)
    to_ret = [WeierstrassUtil._filter_single_landscape(l,bins) 
              for l in landscapes]
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
    
    
def single_area_landscape_bootstrap(area,slice_tmp,skip,N_boostraps):
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
    return heatmap_data,iwt_tmp                                             
    
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
    raw_area_slices = []
    area_bounds = [get_area_bounds(raw_data,area) for area in areas]
    # fix all the spring constants. XXX need to account for this...
    mean_spring_constant = np.mean([r.LowResData.meta.SpringConstant 
                                    for r in raw_data])
    for r in raw_data:
        r.LowResData.meta.SpringConstant = mean_spring_constant
    # fix the constant offset added        
    for i,r in enumerate(raw_data):
        r.Force -= offset_N
    # make the slice needed (the 'full' dataset)
    for i,r in enumerate(raw_data):
        idx_0,idx_f = area_bounds[0][i]
        s = slice(idx_0,idx_f,1)
        this_area = FEC_Util.MakeTimeSepForceFromSlice(r,s)
        raw_area_slices.append(this_area)
        # r is no longer needed; stop referencing it to make space
        raw_data[i] = None
    to_ret = []
    skip = 0
    N_boostraps = 10
    heatmap_data,iwt_tmp = \
        single_area_landscape_bootstrap(areas[0],raw_area_slices,
                                        skip,N_boostraps)   
    filtered_iwt = filter_landscapes(iwt_tmp,area.n_bins)                                        
    return cacheable_data(filtered_iwt,*heatmap_data)                           