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
    for area,slice_tmp in zip(areas,raw_area_slices):
        # get the heatmap histograms
        heatmap_data = get_heatmap_data(slice_tmp,)  
        # get the landscapes 
        iwt_objs = IWT_Util.convert_to_iwt(slice_tmp)
        iwt_helix_data_tmp = \
            InverseWeierstrass.FreeEnergyAtZeroForce(iwt_objs,area.n_bins)
        # make the object we want for this 'area' slice
        to_ret.append(cacheable_data(iwt_helix_data_tmp,heatmap_data))
    return to_ret
    
def make_heatmap(histogram, x_edges,y_edges):
    # XXX ? digitize all the ids so we know what bin they fall into...
    X,Y = np.meshgrid(x_edges,y_edges)
    plt.gca().pcolormesh(X,Y,histogram,cmap=plt.cm.afmhot)
    
def get_energy_landscape_data(data_to_plot,nanometers_to_amino_acids,
                              kT=4.1e-21):
    # get the landscape
    landscape_kT = data_to_plot.landscape.EnergyLandscape/kT
    landscape_kcal_per_mol = landscape_kT * IWT_Util.kT_to_kcal_per_mol()
    extension_nm = data_to_plot.landscape.Extensions*1e9
    # get the change in energy per unit distance (force)
    delta_landscape_kT_per_nm = \
        np.gradient(landscape_kT)/np.gradient(extension_nm)
    delta_landscape_kcal_per_mol_per_nm = \
        delta_landscape_kT_per_nm * IWT_Util.kT_to_kcal_per_mol()
    delta_landscape_kcal_per_mol_per_amino_acid = \
        delta_landscape_kcal_per_mol_per_nm * 1/(nanometers_to_amino_acids)
    return extension_nm,landscape_kcal_per_mol,\
        delta_landscape_kcal_per_mol_per_amino_acid
    
def plot_landscape(extension_nm,landscape_kcal_per_mol,
                   delta_landscape_kcal_per_mol_per_amino_acid,xlim):
    ax_energy = plt.gca()
    plt.plot(extension_nm,landscape_kcal_per_mol,color='k')
    # make a second axis for the number of ammino acids 
    units_energy = r"($\frac{\mathrm{kcal}}{\mathrm{mol}}$)"
    units_energy_delta = r"($\frac{\mathrm{kcal}}{\mathrm{mol} \cdot AA}$)"
    PlotUtilities.lazyLabel("Extension (nm)","Free energy " + units_energy,"")    
    limits_delta = [min(delta_landscape_kcal_per_mol_per_amino_acid),
                    max(delta_landscape_kcal_per_mol_per_amino_acid)]
    plt.xlim(xlim)                    
    label = "Free energy difference " + units_energy_delta
    ax_2 = PlotUtilities.secondAxis(ax_energy,
                                    label=label,color='r',
                                    limits=limits_delta,secondY =True)
    ax_2.plot(extension_nm,delta_landscape_kcal_per_mol_per_amino_acid,
              color='r',linestyle='-',linewidth=0.5)                               
    plt.xlim(xlim)     

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
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    flickering_dir = "../Data/"
    # XXX use the flickering dir for stuff
    cache_dir = flickering_dir 
    GenUtilities.ensureDirExists(flickering_dir)
    meters_to_amino_acids = 1/(0.3e-9)
    nanometers_to_amino_acids = meters_to_amino_acids * 1/1e9
    n_bins = 150
    # write down the areas we want to look at 
    areas = [\
        slice_area([18e-9,75e-9],"Full (no adhesion)",n_bins),
        slice_area([20e-9,27e-9],"Helix A",n_bins),
        slice_area([50e-9,75e-9],"Helix E",n_bins),
        ]    
    # read in the data 
    data_to_analyze = CheckpointUtilities.\
        getCheckpoint("./cached_landscapes.pkl",get_cacheable_data,False,areas,
                      flickering_dir)
    data_to_plot = data_to_analyze[0]         
    heatmap_data = data_to_plot.heatmap_data
    extension_nm,landscape_kcal_per_mol, delta_energy_kcal_per_mol_per_aa = \
            get_energy_landscape_data(data_to_plot,nanometers_to_amino_acids)
    fig = PlotUtilities.figure((3.25,7))    
    # # ploy the heat map 
    ax_heat = plt.subplot(2,1,1)
    heatmap_plot(heatmap_data,nanometers_to_amino_acids)
    xlim_fec = plt.xlim()
    # # plot the energy landscape...
    ax_energy = plt.subplot(2,1,2)    
    plot_landscape(extension_nm,landscape_kcal_per_mol,
                   delta_energy_kcal_per_mol_per_aa,xlim_fec)
    PlotUtilities.savefig(fig,"out.png")

    
if __name__ == "__main__":
    run()
