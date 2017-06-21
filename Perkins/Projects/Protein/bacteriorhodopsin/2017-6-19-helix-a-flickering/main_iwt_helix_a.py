# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re,copy

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import CheckpointUtilities,PlotUtilities,GenUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao
from FitUtil.WormLikeChain.Python.Code import WLC 
import cProfile 

class slice_area:
    def __init__(self,ext_bounds,plot_title,n_bins):
        self.ext_bounds = ext_bounds
        self.plot_title = plot_title
        self.n_bins = n_bins
    @property
    def save_name(self):
        return self.plot_title.replace(" ","_") + ".pkl"

    
def fmt_iwt():
    PlotUtilities.xlabel("Extension (nm)")
    
   
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "./data_in_full/"
    out_dir = "./out_full/"
    GenUtilities.ensureDirExists(out_dir)
    force_read_data = False
    force_iwt = True
    raw_data = IoUtilHao.read_and_cache_data_hao(in_dir,force=force_read_data,
                                                 limit=10)
    adhesion_end_m = 20e-9         
    offset_force_N = 7.1e-12
    min_max_full = [adhesion_end_m,75e-9]
    n_bins_zoom = 100
    n_bins_helix_a = 150
    array = [\
        slice_area([20e-9,27e-9],"Helix A",n_bins_helix_a),
        slice_area([50e-9,75e-9],"Helix E",n_bins_zoom),
        slice_area([57e-9,70e-9],"Helix E (detailed)",n_bins_zoom)
        ]
    # # get the slice we care about for each...
    sliced_data = [ [] for _ in array]
    full_data = []
    for i,d in enumerate(raw_data):
        d.Force -= offset_force_N 
        sliced_fec = FEC_Util.slice_by_separation(d,*min_max_full)
        full_data.append(sliced_fec)
        for i,a in enumerate(array):
            slice_tmp = FEC_Util.slice_by_separation(d,*a.ext_bounds)
            sliced_data[i].append(slice_tmp)
    # # get the IWT of both regions
    n_bins = 100
    n_bins_helix=200
    n_bins_helix_zoom = 60
    iwt_f = InverseWeierstrass.FreeEnergyAtZeroForce
    # convert into the iwt objects needed 
    iwt_full_data= IWT_Util.convert_to_iwt(full_data)
    # get the proper landscapes
    iwt_full = CheckpointUtilities.getCheckpoint(out_dir + "iwt_full.pkl",iwt_f,
                                                 force_iwt,
                                                 iwt_full_data,
                                                 NumBins=n_bins)
    iwt_helices = []
    for i,a in enumerate(array):
        save_name = (out_dir + a.save_name)
        iwt_helix_data_tmp =  IWT_Util.convert_to_iwt(sliced_data[i])
        iwt_tmp = CheckpointUtilities.getCheckpoint(save_name,iwt_f,
                                                    force_iwt,
                                                    iwt_helix_data_tmp,
                                                    NumBins=a.n_bins)  
        iwt_helices.append(iwt_tmp)                                                  
    # make the iwt (energy landscape) plot of the entire protein
    fig = PlotUtilities.figure()
    obj_full_plot = IWT_Plot.plot_free_landscape(iwt_full)    
    plt.axvspan(0,min(obj_full_plot.landscape_ext_nm),alpha=0.3,color='r',
                label="Adhesion region",hatch='/')
    PlotUtilities.title("Full bacteriorhodopsin energy landscape")
    fmt_iwt()
    PlotUtilities.legend()    
    PlotUtilities.savefig(fig,out_dir + "iwt_full.png")
    # plot each of the supregions 
    for i,(a,data) in enumerate(zip(array,sliced_data)):
        fig = PlotUtilities.figure()
        IWT_Plot.plot_free_landscape(iwt_helices[i])    
        fmt_iwt()    
        PlotUtilities.title(a.plot_title)
        PlotUtilities.savefig(fig,out_dir + a.save_name + "_iwt_.png")
    # # make the plots we want
    # make a heat map of all the data and the 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(full_data)
    PlotUtilities.savefig(fig,out_dir + "heat_map.png")
    # plot the zoomed regions 
    for i,(a,data) in enumerate(zip(array,sliced_data)):
        fig = PlotUtilities.figure()
        FEC_Plot.heat_map_fec(sliced_data[i])    
        PlotUtilities.title(a.plot_title)
        PlotUtilities.savefig(fig,out_dir + a.save_name + "_heat_.png")
    # plot each fec individually 
    to_x = lambda x : x*1e9
    to_y = lambda y : y*1e12
    return 
    for i,d in enumerate(raw_data):
        fig = PlotUtilities.figure()
        plt.plot(to_x(d.Separation),to_y(d.Force),color='k',alpha=0.3)
        sep_sliced = full_data[i].Separation
        FEC_Plot._fec_base_plot(to_x(sep_sliced),to_y(full_data[i].Force),
                                style_data=dict(color='r',alpha=0.3))
        plt.xlim([0,2*to_x(max(sep_sliced))])
        plt.ylim([-50,300])
        file_n = GenUtilities.file_name_from_path(d.Meta.SourceFile)
        PlotUtilities.savefig(fig,out_dir + "out{:d}_{:s}.png".format(i,file_n))

if __name__ == "__main__":
    run()
