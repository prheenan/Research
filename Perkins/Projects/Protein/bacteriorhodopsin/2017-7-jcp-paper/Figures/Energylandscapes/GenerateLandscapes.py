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
    def __init__(self,landscape,heatmap_data,heatmap_data_z,k_arr,n_k_arr,
                 originals):
        """
        :param landscape: list, size N, element is list of landscapes with
        spring constant k_arr[i]

        :param heatmap_data</_z>: the heatmap data as a function of extension
        and z (stage position)

        :param k_arr: the spring constant list, size N
        :param n_k_arr:  the number of fecs with k=k_arr[i].
        """
        self.landscape = landscape
        self.heatmap_data = heatmap_data
        self.heatmap_data_z = heatmap_data_z
        self.k_arr = k_arr
        self.n_k_arr = n_k_arr
        self.originals= originals
    def generate_landscape_obj(self):
        """
        Returns: a landscape_data, with all the weights it needs
        """
        all_data = [d for list_v in self.landscape for d in list_v]
        weights = []
        #n_k_arr[i] is the number of fecs used to calculate
        # each landscape in the list of landscapes self.landscape[i]
        data_to_use = []
        original_data_to_use = []
        original_weights = []
        min_n = 20
        for i,list_v in enumerate(self.landscape):
            n_tmp = self.n_k_arr[i]
            if n_tmp < min_n:
                continue
            weights_tmp = [n_tmp for j in range(len(list_v))]
            weights.extend(weights_tmp)
            data_to_use.extend(list_v)
            original_data_to_use.append(self.originals[i])
            original_weights.append(n_tmp)
        # push all the landscapes onto the same grid
        num_bins = [l.q.size for l in data_to_use]
        assert len(set(num_bins)) == 1
        # POST: only one size of landscape
        to_ret = landscape_data(data_to_use,weights=weights,
                                original_data=original_data_to_use,
                                original_weights=original_weights)
        return to_ret
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
    @property
    def dx(self):
        dx = self._extensions_m[1] -  self._extensions_m[0]
        return dx
    def _d_energies_d_m(self,energies):
        return np.array([np.gradient(e) / self.dx for e in energies])
    def _d2_energies_dm2(self,energies):
        return [np.gradient(e)/self.dx 
                for e in (self._d_energies_d_m(energies))]
    def __init__(self,landscape_objs,weights,original_data,original_weights):
        assert len(landscape_objs) > 0 
        self.landscape_objs = landscape_objs
        self.weights = weights
        self.original_weights = original_weights
        self.kT = [1/o.beta for o in landscape_objs][0]
        min_v = min([min(l.q) for l in landscape_objs])
        max_v = max([max(l.q) for l in landscape_objs])
        sizes = np.array([l.q.size for l in landscape_objs])
        n = np.max(sizes)        
        expected_sizes = np.ones(sizes.size)*n
        np.testing.assert_allclose(sizes,expected_sizes)
        self._extensions_m = np.linspace(min_v,max_v,n*10,endpoint=True)
        # save the origial (ie: same sized grid)
        self._extensions_m_original = np.linspace(min_v,max_v,n,endpoint=True)
        self._energies = [l.spline_fit.y(self._extensions_m)
                          for l in landscape_objs]
        self._energies_orig = [l.spline_fit.y(self._extensions_m) for l in
                               original_data]
        # align all the arbitrary zeros
        self._energies = [e - min(e) for e in self._energies]
        self._energies_orig = [e - min(e) for e in self._energies_orig]
        # align the minimums
        self._extensions_m -= min_v
        self._extensions_m_original -= min_v
    def from_Joules_to_kcal_per_mol(self):
        return IWT_Util.kT_to_kcal_per_mol() * (1/self.kT)
    def _raw_uninterpolared_landscapes_kcal_per_mol(self,l):
        return l.G_0 * self.from_Joules_to_kcal_per_mol()
    def amino_acids_per_nm(self):
        return 3     
    def _grid_property(self,f):
        to_ret =[f(l) for l in self.landscape_objs]
        return to_ret 
    @property
    def _extension_grid_nm(self):
        return self._extensions_m * 1e9
    def _get_energy(self,bootstrap):
        return self._energies if bootstrap else self._energies_orig
    def _landscapes_kcal_per_mol(self,bootstrap=True):
        x = self._get_energy(bootstrap)
        return np.array(x) * self.from_Joules_to_kcal_per_mol()
    def _delta_landscapes_kcal_per_mol_per_nm(self,bootstrap=True):
        x = self._get_energy(bootstrap)
        energy_kcal_per_mol_per_m = \
            self._d_energies_d_m(x) * self.from_Joules_to_kcal_per_mol()
        energy_kcal_per_mol_per_nm = energy_kcal_per_mol_per_m * 1e-9
        return energy_kcal_per_mol_per_nm
    def _delta_landscapes_kcal_per_mol_p_AA(self,**kw):
        energy_kcal_per_mol_per_AA = \
            self._delta_landscapes_kcal_per_mol_per_nm(**kw)/self.amino_acids_per_nm()
        return energy_kcal_per_mol_per_AA
    def mean_std_opt(self,weights):
        return dict(axis=0,weights=weights)
    def _avg(self,x,**kw):
        return np.average(x,**self.mean_std_opt(**kw))
    def _std(self,x,weights):
        mean_tmp = self._avg(x,weights=weights)
        variance = self._avg((x-mean_tmp)**2,weights=weights)
        return np.sqrt(variance)
    @property        
    def mean_landscape_kcal_per_mol(self):
        # get the landscapes, XXX need to interpolate back onto uniform grid.
        to_avg = self._landscapes_kcal_per_mol(bootstrap=False)
        return self._avg(to_avg,weights=self.original_weights)
    @property                
    def std_landscape_kcal_per_mol(self):
        return self._std(self._landscapes_kcal_per_mol(bootstrap=True),
                         weights=self.weights)
    @property        
    def mean_delta_landscape_kcal_per_mol_per_nm(self):
        to_avg = self._delta_landscapes_kcal_per_mol_per_nm(bootstrap=False)
        return self._avg(to_avg,weights=self.original_weights)
    @property                               
    def std_delta_landscape_kcal_per_mol_per_nm(self):
        to_std = self._delta_landscapes_kcal_per_mol_per_nm(bootstrap=True)
        return self._std(to_std,weights=self.weights)

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
    
    
def filter_landscapes(landscapes,n_bins,bins=None,**kw):
    # zero everything; don't care about absolute offsets.
    min_v = max([min(l.q) for l in landscapes])
    max_v = min([max(l.q) for l in landscapes])
    if (bins is None):
        bins = np.linspace(min_v,max_v,n_bins)
    to_ret = [WeierstrassUtil._filter_single_landscape(l,bins,**kw)
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
    
def get_heatmap_data(time_sep_force_arr,bins=(300,100)):
    sep_nm = [t.Separation*1e9 for t in time_sep_force_arr]
    z_nm = [t.ZSnsr*1e9 for t in time_sep_force_arr]
    force_pN = [t.Force*1e12 for t in time_sep_force_arr]
    heatmap_force_extension = heatmap(sep_nm,force_pN,bins=bins)
    heatmap_force_z = heatmap(z_nm,force_pN,bins=bins)
    return heatmap_force_extension,heatmap_force_z
    
    
def single_area_landscape_bootstrap(area,slice_tmp,skip,N_boostraps,cache_dir):
    # get the landscapes (we use N, to get an error)
    iwt_objs = WeierstrassUtil.convert_list_to_iwt(slice_tmp)
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
    GenUtilities.ensureDirExists(cache_dir)
    iwt_tmp = CheckpointUtilities.multi_load(cache_dir,load_func=functor,
                                             force=False,
                                             name_func=name_func)   
    return iwt_tmp
    
def free_list(l):
    for l_tmp in l:
        l_tmp = None
    del l
    
def get_raw_data(flickering_dir,areas,heat_bins=(100,100),
                  offset_N=7.1e-12):
    raw_data = IoUtilHao.read_and_cache_data_hao(None,force=False,
                                                 cache_directory=flickering_dir,
                                                 limit=100,
                                                 renormalize=False)
    # only look at data with ~300nm/s
    v_exp = 300e-9
    raw_data = [r for r in raw_data 
                if np.allclose(r.Velocity,v_exp,atol=0,rtol=0.01)]
    raw_area_slices = []
    area_bounds = [get_area_bounds(raw_data,a) for a in areas]
    # fix the constant offset added, see SI of Science paper
    for i,r in enumerate(raw_data):
        r.Force -= offset_N
    # make the slice needed (the 'full' dataset)
    raw_area_slice = []
    for i,r in enumerate(raw_data):
        idx_0,idx_f = area_bounds[0][i]
        s = slice(idx_0,idx_f,1)
        this_area = FEC_Util.MakeTimeSepForceFromSlice(r,s)
        raw_area_slice.append(this_area)
        # r is no longer needed; stop referencing it to make space
        raw_data[i] = None
    # get the heatmap on the entire slice
    heatmap_data = get_heatmap_data(raw_area_slice)
    min_data = 20
    area_of_interest = areas[0]
    k_arr_raw = [r.LowResData.meta.SpringConstant for r in raw_area_slice]
    k_set = np.array(sorted(list(set(k_arr_raw))))
    k_idx = np.array([np.argmin(np.abs(k_tmp - k_set)) for k_tmp in k_arr_raw])
    data_to_use = []
    k_arr = []
    for i in range(len(k_idx)):
        tmp_idx = np.where(np.abs(k_idx - i) < 1e-6)[0]
        m_data = [raw_area_slice[j] for j in tmp_idx]
        len_data = len(m_data)
        # record the actual spring constants and data to use. 
        if (len_data >= min_data):
            data_to_use.append(m_data)
            k_arr.append(k_set[i])
    for d,k_tmp in zip(data_to_use,k_arr):
        yield RawData(d,k_tmp,k_set,heatmap_data,k_arr)
        
class RawData:
    def __init__(self,data,k_tmp,k_set,heatmap,k_arr):
        self.data=data
        self.k_tmp = k_tmp
        self.k_set = k_set
        self.heatmap = heatmap
        self.k_arr = k_arr
    
def get_filtered_original_landscape(fecs,n_bins):
    iwt_tmp_orig =  WeierstrassUtil.convert_list_to_iwt(fecs)
    # get the landscape and filter it.
    original = \
        InverseWeierstrass.free_energy_inverse_weierstrass(iwt_tmp_orig)     
    original_filtered = \
        filter_landscapes([original],n_bins)
    # free the original data
    free_list(iwt_tmp_orig)        
    # append the first (only) element
    assert len(original_filtered) == 1        
    return original_filtered[0]     
    
def get_cacheable_data(areas,flickering_dir,cache_dir="./cache"):
    # save out the raw data we need 
    area_of_interest = areas[0]    
    functor = lambda : get_raw_data(flickering_dir,areas)
    name_func = lambda  i_tmp,d: \
        area_of_interest.save_name + "_raw_{:d}".format(i_tmp) 
    # get a generator (!) for all the raw data...
    data = CheckpointUtilities.multi_load(cache_dir="./raw_data/",
                                          load_func=functor,
                                          force=False,
                                          name_func=name_func)      
    # POST: data_to_use[i] has the data for the spring constant i
    # POST: all data are sorted by spring constant
    data_lengths = []
    filtered_iwt = []
    originals = []
    N_boostraps = 10
    skip = 0    
    for i,d in enumerate(data):
        data_lengths.append(len(d.data))
        original_name = "./original_{:.3g}pN_nm.pkl".format(1000*d.k_tmp)
        original_filtered = \
            CheckpointUtilities.getCheckpoint(original_name,
                                              get_filtered_original_landscape,
                                              False,
                                              d.data,area_of_interest.n_bins)
        cache_dir = "./{:s}_k_{:.4g}pN_nm/".\
                format(area_of_interest.save_name,d.k_set[i]*1000)
        iwt_tmp = single_area_landscape_bootstrap(area_of_interest,d.data,
                                                  skip,N_boostraps,
                                                  cache_dir=cache_dir)
        # immediately filter; don't use the unfiltered data. avoids memory
        # problems.
        filtered_tmp = filter_landscapes(iwt_tmp,area_of_interest.n_bins)
        filtered_iwt.append(filtered_tmp)
        # explicitly set iwt_tmp to None, to have it be garbage collected
        # get the landscape using all the original data
        # get the filtered version...
        free_list(iwt_tmp)
        originals.append(original_filtered)
    # POST: filtered_iwt[i] has all bootstraps from fecs with spring k_set[i]
    key_d = d
    assert len(data_lengths) > 0 , "Couldn't find any data to fit"    
    return cacheable_data(filtered_iwt,*key_d.heatmap,k_arr=key_d.k_arr,
                          n_k_arr=data_lengths,originals=originals)
