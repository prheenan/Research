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
    
def get_cacheable_data(areas,flickering_dir):
    force_read_data = False    
    raw_data = IoUtilHao.read_and_cache_data_hao(None,force=force_read_data,
                                                 cache_directory=flickering_dir,
                                                 limit=None)
    raw_area_slices = []
    for area in areas:
        this_area = [FEC_Util.slice_by_separation(r,*area.ext_bounds) 
                     for r in raw_data]
        raw_area_slices.append(this_area)
    to_ret = []
    for slice_tmp in raw_area_slices:
        heatmap_data = get_heatmap_data(slice_tmp,bins=(100,100))  
        to_ret.append(cacheable_data(None,heatmap_data))
    return to_ret
    
def make_heatmap(histogram, x_edges,y_edges):
    # XXX ? digitize all the ids so we know what bin they fall into...
    X,Y = np.meshgrid(x_edges,y_edges)
    plt.gca().pcolormesh(X,Y,histogram,cmap=plt.cm.afmhot)
    
    
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
    meters_to_amino_acids = 1/0.3
    n_bins = 150
    areas = [\
        slice_area([18e-9,75e-9],"Full (no adhesion)",n_bins),
        slice_area([20e-9,27e-9],"Helix A",n_bins),
        slice_area([50e-9,75e-9],"Helix E",n_bins),
        slice_area([57e-9,70e-9],"Helix E (detailed)",n_bins)
        ]    
    data_to_analyze = CheckpointUtilities.\
        getCheckpoint("./cached_landscapes.pkl",get_cacheable_data,False,areas,
                      flickering_dir)
    heatmap_data = data_to_analyze[0].heatmap_data
    fig = PlotUtilities.figure((3.25,7))
    ax_heat = plt.subplot(2,1,1)
    make_heatmap(*heatmap_data)
    PlotUtilities.lazyLabel("","Force (pN)","")
    # make a second axis for the number of ammino acids 
    xlim_fec = plt.xlim()
    limits = np.array(xlim_fec) * meters_to_amino_acids
    PlotUtilities.secondAxis(ax_heat,"Extension (AA #)",limits,secondY =False)
    plt.xlim(xlim_fec)
    energy_axis = plt.subplot(2,1,2)    
    PlotUtilities.lazyLabel("","Free energy (kT)","")
    # make a second axis for the number of ammino acids 
    IWT_Plot.format_kcal_per_mol_second_axis_after_kT_axis(ax=energy_axis)
    plt.xlim(xlim_fec)
    PlotUtilities.savefig(fig,"out.png")

    
if __name__ == "__main__":
    run()
