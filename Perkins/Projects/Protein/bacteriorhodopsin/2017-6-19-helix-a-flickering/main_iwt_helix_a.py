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
   
def convert_to_iwt(time_sep_force,frac_vel=0.1):
    iwt_data = [IWT_Util.ToIWTObject(d) for d in time_sep_force]
    set_vel = IWT_Util.set_separation_velocity_by_first_frac    
    for d in iwt_data:
        set_vel(d,fraction_for_vel=frac_vel)  
    return iwt_data
    
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
    in_dir = "./data_in/"
    out_dir = "./out/"
    GenUtilities.ensureDirExists(out_dir)
    force_read_data = False
    force_iwt = False
    data = IoUtilHao.read_and_cache_data_hao(in_dir,force=force_read_data,
                                             limit=40)
    adhesion_end_m = 20e-9                                             
    min_max_full = [adhesion_end_m,75e-9]
    min_max_helix_a = [45e-9,75e-9]
    min_max_helix_a_zoomed = [57e-9,61e-9]
    # # get the slice we care about for each...
    sliced_data = []
    helix_a_data = []
    helix_a_zoomed_data = []
    for i,d in enumerate(data):
        sliced_fec = FEC_Util.slice_by_separation(d,*min_max_full)
        slice_helix_a = FEC_Util.slice_by_separation(d,*min_max_helix_a)
        slice_helix_a_zoomed = \
            FEC_Util.slice_by_separation(d,*min_max_helix_a_zoomed)
        sliced_data.append(sliced_fec)
        helix_a_data.append(slice_helix_a)
        helix_a_zoomed_data.append(slice_helix_a_zoomed)
    # # get the IWT of both regions
    n_bins = 100
    n_bins_helix=200
    n_bins_helix_zoom = 60
    iwt_f = InverseWeierstrass.FreeEnergyAtZeroForce
    # convert into the iwt objects needed 
    iwt_full_data= convert_to_iwt(sliced_data)
    iwt_helix_a_data =  convert_to_iwt(helix_a_data)
    iwt_helix_a_zoomed_data = convert_to_iwt(helix_a_zoomed_data)
    # get the proper landscapes
    iwt_full = CheckpointUtilities.getCheckpoint("iwt_full.pkl",iwt_f,force_iwt,
                                                 iwt_full_data,
                                                 NumBins=n_bins)
    iwt_helix_a = CheckpointUtilities.getCheckpoint("iwt_helix_a.pkl",iwt_f,
                                                    force_iwt,
                                                    iwt_helix_a_data,
                                                    NumBins=n_bins_helix)  
    iwt_helix_a_zoomed = \
        CheckpointUtilities.getCheckpoint("iwt_helix_a_zoomed.pkl",iwt_f,
                                          force_iwt,
                                          iwt_helix_a_zoomed_data,
                                          NumBins=n_bins_helix_zoom)                                                      
    # make the iwt (energy landscape) plot of the entire protein
    fig = PlotUtilities.figure()
    obj_full_plot = IWT_Plot.plot_free_landscape(iwt_full)    
    plt.axvspan(0,min(obj_full_plot.landscape_ext_nm),alpha=0.3,color='r',
                label="Adhesion region",hatch='/')
    PlotUtilities.title("Full bacteriorhodopsin energy landscape")
    fmt_iwt()
    PlotUtilities.legend()    
    PlotUtilities.savefig(fig,"./iwt_full.png")
    # make the iwt (energy landscape) plot of the helix, including a     
    fig = PlotUtilities.figure()
    IWT_Plot.plot_free_landscape(iwt_helix_a)    
    fmt_iwt()    
    PlotUtilities.savefig(fig,"./iwt_helix_a.png") 
    # make the zoomed plot of the helix a region 
    fig = PlotUtilities.figure()
    IWT_Plot.plot_free_landscape(iwt_helix_a_zoomed)    
    fmt_iwt()    
    PlotUtilities.savefig(fig,"./iwt_helix_a_zoomed.png")     
    # # make the plots we want
    # make a heat map of all the data and 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(sliced_data)
    PlotUtilities.savefig(fig,out_dir + "heat_map.png")
    # make a heat map of just the helix a region
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(helix_a_data)
    PlotUtilities.savefig(fig,out_dir + "heat_map_helix_a.png")    
    # and of the zoomed region 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(helix_a_zoomed_data)
    PlotUtilities.savefig(fig,out_dir + "heat_map_helix_a_zoom.png")    
    # plot each individually 
    to_x = lambda x : x*1e9
    to_y = lambda y : y*1e12
    for i,d in enumerate(data):
        fig = PlotUtilities.figure()
        plt.plot(to_x(d.Separation),to_y(d.Force),color='k',alpha=0.3)
        sep_sliced = sliced_data[i].Separation
        FEC_Plot._fec_base_plot(to_x(sep_sliced),to_y(sliced_data[i].Force),
                                style_data=dict(color='r',alpha=0.3))
        plt.xlim([0,2*to_x(max(sep_sliced))])
        plt.ylim([-50,300])
        file_n = GenUtilities.file_name_from_path(d.Meta.SourceFile)
        PlotUtilities.savefig(fig,out_dir + "out{:d}_{:s}.png".format(i,file_n))

if __name__ == "__main__":
    run()
