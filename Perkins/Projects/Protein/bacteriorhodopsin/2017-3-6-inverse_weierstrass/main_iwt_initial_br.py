# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re,copy

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from Research.Perkins.Projects.Protein.bacteriorhodopsin import IoUtilHao
from FitUtil.WormLikeChain.Python.Code import WLC 
import cProfile 
   
def get_first_peak_slices(absolute_data_dir):
    downsample_n = 1
    fraction_for_vel = 0.1    
    limit = 20
    force_sample = False
    force= False
    cache_directory = "./"
    out_base = cache_directory    
    retracts = IoUtilHao.get_downsampled_data(downsample_n,force_sample,
                                              cache_directory,
                                              absolute_data_dir,limit=limit)
    data_fec = [IoUtilHao.get_retract_pulling_region(d,zero=True) 
                for d in retracts]
    # get a x% 'fudge factor' used for filtering etc.
    frac = 0.02
    fudge_N = [int(np.ceil(x.Force.size*frac)) for x in data_fec]
    filtered = [FEC_Util.GetFilteredForce(x,n)
                for x,n in zip(data_fec,fudge_N)]
    slice_max = lambda x: slice(0,np.argmax(x.Force),1)
    # get the slices until the max (used to fit to the FEC)
    max_slices = [slice_max(r_filtered) for r_filtered in filtered]
    until_max = [FEC_Util.MakeTimeSepForceFromSlice(r,slice(0,slice_v.stop))
                 for r,slice_v in zip(data_fec,max_slices)]
    for i,(data_slice,full) in enumerate(zip(until_max,data_fec)):         
        fig = PlotUtilities.figure()
        plt.plot(data_slice.Separation*1e9,data_slice.Force*1e12,color='r')
        plt.plot(full.Separation*1e9,full.Force*1e12,color='k',alpha=0.3)
        PlotUtilities.lazyLabel("Sep (nm)","Force (pN)","")
        PlotUtilities.savefig(fig,"./out/full_fec_{:d}.png".format(i))
    # convert all the FEC into IWT                                             
    data_iwt = []
    for r in until_max:
        tmp_iwt = IWT_Util.ToIWTObject(r)
        # set all the effective velocities
        IWT_Util.set_separation_velocity_by_first_frac(tmp_iwt,fraction_for_vel)
        data_iwt.append(tmp_iwt)     
    return data_iwt
    
def get_contour_lengths(data_iwt,kwargs_fit):
    contour_lengths = []
    for reference in data_iwt:
        max_ext = max(reference.Separation)
        brute_dict = dict(ranges=[ [max_ext/2,max_ext*2] ],Ns=15)
        sep,force = reference.Separation,reference.Force                      
        x0,y = WLC.fit(sep,force,brute_dict=brute_dict,
                       **kwargs_fit)
        contour_lengths.append(x0)
    return np.concatenate(contour_lengths)
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = FEC_Util.default_data_root()
    # XXX use Haos...
    absolute_data_dir = "./data_in/"
    out_base = "./"
    fraction_for_vel = 0.1
    # get the region near the first peak
    data_iwt = \
        CheckpointUtilities.getCheckpoint("./peaks.pkl",
                                          get_first_peak_slices,
                                          False,
                                          absolute_data_dir)                                   
    kwargs_fit = dict(kbT = 4.11e-21,
                      Lp=0.3e-9,
                      K0=1000e-12)                            
    # get the contour lengths                                                 
    contour_lengths = \
        CheckpointUtilities.getCheckpoint("./contour.pkl",
                                          get_contour_lengths,False,
                                          data_iwt,
                                          kwargs_fit)                                            
    # align all of them so that the contour length is at 0+x
    offset_nm = 0
    iwt_offset = []
    for d,L0 in zip(data_iwt,contour_lengths):
        tmp = copy.deepcopy(d)
        tmp.Extension -= (L0 - offset_nm)
        iwt_offset.append(tmp)
    # plot *all* of the data, along with the WLC fits...
    for i,(r,x0) in enumerate(zip(data_iwt,contour_lengths)):
        sep,force = r.Separation,r.Force
        grid_x,grid_y,predicted_force = \
            WLC.inverted_wlc(sep,force,x0,**kwargs_fit)   
        fig = PlotUtilities.figure()
        FEC_Plot._fec_base_plot(sep*1e9,force*1e12)
        plt.plot(grid_x*1e9,grid_y*1e12,color='r')
        PlotUtilities.lazyLabel("Separation","Force","")
        PlotUtilities.savefig(fig,"./{:d}_FEC.png".format(i))
    # plot the data before aligning
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(data_iwt)
    PlotUtilities.savefig(fig,"heatmap_before.png")                 
    # plot the data after aligning
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(iwt_offset)
    PlotUtilities.savefig(fig,"heatmap_after.png")       
    # POST: they are all set. get the IWT 
    num_bins = 250
    pr = cProfile.Profile()
    pr.enable()
    LandscapeObj =  InverseWeierstrass.\
        FreeEnergyAtZeroForce(iwt_offset,NumBins=num_bins)
    pr.disable()
    pr.print_stats(sort='time') 
    force_N = [1e-12,5e-12,10e-12,20e-12,40e-12,100e-12,250e-12,500e-12]
    max_extension_nm = max([max(f.Separation) for f in iwt_offset]) * 1e9
    ylim_kT = [-5,None]
    for i,f in enumerate(force_N):  
        xlim_nm = [min(LandscapeObj.Extensions)*1e9,
                   max_extension_nm]
        fig = PlotUtilities.figure(figsize=(6,12))
        kw_landscape = dict(f_one_half_N=f)
        # # plot the heat map after aligning
        plt.subplot(3,1,1)
        FEC_Plot.heat_map_fec(iwt_offset,use_colorbar=False)
        PlotUtilities.xlabel("")
        plt.xlim(xlim_nm)
        # # plot the free and tilted landscape
        ax_free = plt.subplot(3,1,2)
        IWT_Plot.plot_free_landscape(LandscapeObj,**kw_landscape)  
        PlotUtilities.xlabel("")
        plt.xlim(xlim_nm)        
        plt.ylim(ylim_kT)
        # add a second axis for kcal/mol
        ylim_kT = np.array(plt.ylim())
        ylim_kcal_mol = IWT_Util.kT_to_kcal_per_mol(ylim_kT)
        PlotUtilities.secondAxis(ax_free,label="kcal/mol",limits=ylim_kcal_mol,
                                 secondY =True)
        ax_tilt = plt.subplot(3,1,3)
        IWT_Plot.plot_tilted_landscape(LandscapeObj,fmt_f_label="{:.0f}",
                                       **kw_landscape) 
        PlotUtilities.xlabel("")
        plt.xlim(xlim_nm)       
        plt.ylim(ylim_kT)        
        # add another second axis for kcal/mol
        PlotUtilities.secondAxis(ax_tilt,label="kcal/mol",limits=ylim_kcal_mol,
                                 secondY =True)
        out_name= out_base + "IWT{:d}_{:.1g}.png".format(i,f*1e12)   
        PlotUtilities.savefig(fig,out_name)


if __name__ == "__main__":
    run()
