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
    ylim_kT = np.array(plt.ylim())
    ylim_kcal_per_mol = IWT_Util.kT_to_kcal_per_mol() * ylim_kT
    PlotUtilities.xlabel("Extension (nm)")
    PlotUtilities.secondAxis(ax=plt.gca(),label="Energy (kcal/mol)",
                             limits=ylim_kcal_per_mol,color='b',secondY=True)
   
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
    data = IoUtilHao.read_and_cache_data_hao(in_dir,force=False,limit=40)
    min_max_no_adhesion = [15e-9,75e-9]
    min_max_helix_a = [55e-9,65e-9]
    # # get the slice we care about for each...
    sliced_data = []
    helix_a_data = []
    for i,d in enumerate(data):
        sliced_fec = FEC_Util.slice_by_separation(d,*min_max_no_adhesion)
        slice_helix_a = FEC_Util.slice_by_separation(d,*min_max_helix_a)
        sliced_data.append(sliced_fec)
        helix_a_data.append(slice_helix_a)
    # # get the IWT of both regions
    n_bins = 100
    n_bins_helix=200
    iwt_f = InverseWeierstrass.FreeEnergyAtZeroForce
    # convert into the iwt objects needed 
    iwt_full_data= convert_to_iwt(sliced_data)
    iwt_helix_a_data =  convert_to_iwt(helix_a_data)
    # get the proper landscapes
    iwt_full = CheckpointUtilities.getCheckpoint("iwt_full.pkl",iwt_f,False,
                                                 iwt_full_data,
                                                 NumBins=n_bins)
    iwt_helix_a = CheckpointUtilities.getCheckpoint("iwt_helix_a.pkl",iwt_f,
                                                    False,
                                                    iwt_helix_a_data,
                                                    NumBins=n_bins_helix)  
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
    fig = PlotUtilities.figure(figsize=(4.5,6))
    IWT_Plot.plot_free_landscape(iwt_helix_a)    
    fmt_iwt()    
    PlotUtilities.xlabel("")
    PlotUtilities.savefig(fig,"./iwt_helix_a.png")  
    # # make the plots we want
    # make a heat map of all the data and 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(sliced_data)
    PlotUtilities.savefig(fig,out_dir + "heat_map.png")
    # make a heat map of just the helix a region
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(helix_a_data)
    PlotUtilities.savefig(fig,out_dir + "heat_map_helix_a.png")    
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
