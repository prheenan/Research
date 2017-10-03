# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,cProfile,os,copy

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from FitUtil.FreelyJointedChain.Python.Code import FJC
from GeneralUtil.python import PlotUtilities,CheckpointUtilities,GenUtilities
from Research.Personal.EventDetection.Util import Analysis
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass,WeierstrassUtil

def hairpin_plots(example,filter_fraction,out_path):
    n_filter = int(np.ceil(example.Force.size * filter_fraction))
    # make a plot vs time 
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    FEC_Plot.force_versus_time(example,NFilterPoints=n_filter)
    plt.subplot(2,1,2)
    FEC_Plot.z_sensor_vs_time(example,NFilterPoints=n_filter)
    PlotUtilities.savefig(fig,out_path + "vs_time.png")
    # make a force-extension plot
    fig = PlotUtilities.figure()
    FEC_Plot.FEC(example,NFilterPoints=n_filter)
    PlotUtilities.savefig(fig,out_path + "vs_sep.png")

def make_energy_landscape_plots(out_dir,energy_landscape_unfolding):
    limits_kT = [-2,15]
    landscape_zeroed = copy.deepcopy(energy_landscape_unfolding)
    landscape_zeroed.Extensions -= min(landscape_zeroed.Extensions)
    # plot the iwt transform as just a free landscape                      
    fig = PlotUtilities.figure()
    IWT_Plot.plot_free_landscape(landscape_zeroed)
    PlotUtilities.xlabel("Separation (nm)")
    PlotUtilities.savefig(fig,out_dir + "free_landscape.png")
    # plot the tileted landscape 
    f_one_half_N_arr = np.array([7,8.5,10,12,13,15])*1e-12        
    for i,tilt_N in enumerate(f_one_half_N_arr):
        fig = PlotUtilities.figure()    
        IWT_Plot.plot_tilted_landscape(landscape_zeroed,f_one_half_N=tilt_N)
        PlotUtilities.xlabel("Separation (nm)")    
        save_name = "{:s}free_landscape_tilted_{:d}_{:.2g}.png".\
                    format(out_dir,i,tilt_N*1e12)
        PlotUtilities.savefig(fig,save_name)
    
    
def fit_polymer_model(example):
    """
    Fits the polymer model de jour to example
    
    Args:
        example: the timesepforce to fit to
    Returns:
        x0, the parameters of the fit 
    """
    wlc_params = dict(K0=2000e-12,kbT=4.1e-21)
    ranges = [(10e-9,90e-9),(0.1e-9,1e-9)]
    fit_dict = dict(brute_dict=dict(Ns=10,ranges=ranges),
                    **wlc_params)
    x_raw,y_raw = example.Separation,example.Force
    x0,model_x,model_y = FJC.fit_fjc_contour(x_raw,y_raw,**fit_dict)        
    return x0
    
def get_polymer_coefficients(split_fecs,working_distance_nm):
    """
    returns the polymer coefficients for the given examples, as well as the 
    'valid' indices (ie: if we couldn't fit fot some reason, that idx wont be 
    there)
    
    Args;
        examples:
        working_distance_nm: currently, we assume     
    Returns:
        tuple of <coefficients, valid indices>
    """
    coeffs,idx_valid = [],[]
    for i,s in enumerate(split_fecs):
        retract = s.retract 
        max_fit_idx = np.argmax(retract.Force)
        # find the first idx we are <= the working distance
        target_sep = (retract.Separation[max_fit_idx] - \
                      working_distance_nm * 1e-9)
        idx = np.where(retract.Separation <= target_sep)[0]
        if (idx.size == 0):
            continue
        # POST: at least one thing to fit 
        fit_slice = slice(idx[-1],max_fit_idx,1)
        retract_fit = FEC_Util.MakeTimeSepForceFromSlice(retract,fit_slice)
        x0 = fit_polymer_model(retract_fit)
        # offset the retract by the contour length (first paramter)
        coeffs.append(x0)
        idx_valid.append(i)
    return coeffs,idx_valid

def relative_idx_start_of_final_retract(split_fec):
    """
    Returns: the index into split_fec.retract where the 'final' approach
             starts (ie: no unfolding/refolding). (Indexes must be in the Note)
    """
    retr= split_fec.retract
    idx = [int(i) for i in retr.Meta.Indexes.split(",")]
    # last index is end, second to last is end of 'dwell' (indenter )
    end_of_unfolding_idx = idx[-2] 
    # subtract off the end of the actual approach and surface dwell
    n_points_approach_dwell = split_fec.n_points_approach_dwell()
    # this is to get the *relative* slice index into retr
    # since retr has index 0 and point n_points_approach_dwell in the fec
    # this *assumes* that the DwellTime is set to the indenter time 
    relative_unfolding_start_idx = end_of_unfolding_idx-n_points_approach_dwell
    return relative_unfolding_start_idx
    
def get_unfolding_slice_only(split_fec):
    """
    See: get_unfolding_and_refolding_slices, except only returns the (final)
    retract, after the indenter is done (ie: no unfolding/refolding stuff)
    
    Args:
        split_fec: see get_unfolding_and_refolding_slices
    returns:
        TimeSepForce with just the unfolding portion.
    """
    retr = split_fec.retract 
    relative_unfolding_start_idx = \
        relative_idx_start_of_final_retract(split_fec)
    slice_unfolding_only = slice(relative_unfolding_start_idx,None,1)
    # slice and return the object 
    unfolding_only = FEC_Util.MakeTimeSepForceFromSlice(retr,
                                                        slice_unfolding_only)
    return unfolding_only
    
def get_unfolding_and_refolding_slice(split_fec):
    """
    Given a split_fec with a single unnfolding/folding curve, returns the 
    approach and retract points, trying to keep the minimum and maximum seps
    consistent. 
    
    Args:
        split_fec: the split force extension curve to use; retract should
        *include* the indenter region and "Indexes" meta properly set 

    Returns:
        just the single approach retract 
    """
    relative_unfolding_start_idx = \
        relative_idx_start_of_final_retract(split_fec)
    retr = split_fec.retract
    slice_refolding_experiment = slice(0,relative_unfolding_start_idx,1)
    # get the region up until the start of the final retract 
    slice_to_use = FEC_Util.\
        MakeTimeSepForceFromSlice(retr,slice_refolding_experiment)
    return slice_to_use
    
        
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    abs_dir = "./"
    cache_dir = "./cache/"
    out_dir = "./out/"
    force_iwt = True 
    GenUtilities.ensureDirExists(cache_dir)
    GenUtilities.ensureDirExists(out_dir)    
    examples = FEC_Util.\
        cache_individual_waves_in_directory(pxp_dir=abs_dir,force=False,
                                            cache_dir=cache_dir,limit=20)
    # filter all the fecs 
    good_splits_original = copy.deepcopy(examples)
    n_points_f = lambda x: int(np.ceil(2e-3*s.Force.size))
    examples = [FEC_Util.GetFilteredForce(s,n_points_f(s))
                for s in examples]                                            
    ### XXX TODO
    # (1) correct for interference artifact
    # (2) get regions for WLC fit
    # (3) fit WLC to regions
    # (4) Invert WLC, determine dsDNA and ssDNA contour lengths at each force 
    # split the fecs...
    split_fecs = []
    for i,ex in enumerate(examples):
        split_fec = Analysis.zero_and_split_force_extension_curve(ex)
        retract = split_fec.retract
        split_fecs.append(split_fec)  
    # align them by the contour lengths (assuming we have at least xnm to 
    # work with 
    working_distance_nm = 30
    coeffs,idx = CheckpointUtilities.\
        getCheckpoint("./polymer.pkl",get_polymer_coefficients,False,
                      split_fecs,working_distance_nm)
    # get the split fecs we could actually fit
    good_splits = [split_fecs[i] for i in idx]
    # align all the retracts by the contour lenghts
    contour_L0 = [c[0] for c in coeffs]
    arbitrary_offset = 90e-9
    for L0,split_fec in zip(contour_L0,good_splits):
        split_fec.retract.Separation -= L0
        split_fec.retract.Separation += arbitrary_offset
    # POST: retracts are all aligned.         
    # for a simple IWT, only look at until the unfolding region
    unfolding_retracts = [get_unfolding_slice_only(split_fec) 
                          for e in good_splits]
    refolding_experiments = \
        [get_unfolding_and_refolding_slice(r) for r in good_splits]
    # get the extension maximum and minimum bounds. 
    ext_min_m = lambda s: min(s.Separation)
    ext_max_m = lambda s: min(s.Separation) + 50e-9
    # slice the refolding experiments 
    sliced_refolds = [FEC_Util.slice_by_separation(s,ext_min_m(s),ext_max_m(s))
                      for s in refolding_experiments]                      
    # split the refolding experiments into iwt       
    iwt_refolds = [ \
        WeierstrassUtil.split_into_iwt_objects(s,
                                        fraction_for_vel=0.1,
                                        f_split=IWT_Util.split_by_max_sep)
        for s in  sliced_refolds]                                                                   
    unfolding_objs = [u[0] for u in iwt_refolds]
    refolding_objs = [u[1] for u in iwt_refolds]   
    # slice to just the first L0 (before the final rupture)
    max_meters = arbitrary_offset
    final_rupture_only = [FEC_Util.slice_by_separation(u,-np.inf,max_meters) 
                          for u in unfolding_retracts]
    # convert to the type iwt needs                          
    final_unfolding_iwt = \
        [WeierstrassUtil.convert_to_iwt(r,frac_vel=0.2) 
         for r in final_rupture_only]
    # get the iwt tx 
    n_bins = 200
    energy_landscape_unfolding = CheckpointUtilities.\
        getCheckpoint("./landscape.pkl",
                      InverseWeierstrass.FreeEnergyAtZeroForce,force_iwt,
                      final_unfolding_iwt,NumBins=n_bins)
    # make a heat map of all retracts 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec([r.retract for r in good_splits])
    PlotUtilities.title("FEC Heat map, aligned by L0, N={:d}".\
                        format(len(good_splits)))
    PlotUtilities.savefig(fig,out_dir + "heat.png")    
    # make a heat map of just the region for the unfolding iwt 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(final_rupture_only)
    PlotUtilities.title("FEC Final unfolding heat map, aligned by L0, N={:d}".\
                        format(len(good_splits)))
    PlotUtilities.savefig(fig,out_dir + "heat_unfolding.png")                           
    # make a heat map of the unfolding and refolding  experiment data 
    fig = PlotUtilities.figure(figsize=(4,7))
    plt.subplot(2,1,1)
    FEC_Plot.heat_map_fec(unfolding_objs)
    plt.subplot(2,1,2)
    FEC_Plot.heat_map_fec(refolding_objs)
    PlotUtilities.savefig(fig,out_dir + "heat_refolding.png")         
    # get the refolding one                      
    energy_landscape_bidirectional_folding = CheckpointUtilities.\
        getCheckpoint("./landscape_bidirectional.pkl",
                      InverseWeierstrass.FreeEnergyAtZeroForce,force_iwt,
                      unfolding_objs,RefoldingObjs=refolding_objs,
                      NumBins=100)                      
    energy_landscape_bidirectional_only_unfolding = CheckpointUtilities.\
        getCheckpoint("./landscape_bidirectional_only_unfold.pkl",
                      InverseWeierstrass.FreeEnergyAtZeroForce,force_iwt,
                      unfolding_objs,RefoldingObjs=[],
                      NumBins=100)   
    energy_landscape_bidirectional_only_refolding = CheckpointUtilities.\
        getCheckpoint("./landscape_bidirectional_only_refold.pkl",
                      InverseWeierstrass.FreeEnergyAtZeroForce,force_iwt,
                      refolding_objs,RefoldingObjs=[],
                      NumBins=100)                           
    make_energy_landscape_plots(out_dir +"bi_",
                                energy_landscape_bidirectional_folding)
    make_energy_landscape_plots(out_dir +"same_",
                                energy_landscape_unfolding)                                
    make_energy_landscape_plots(out_dir +"bi_only_unfold",
                                energy_landscape_bidirectional_only_unfolding)  
    make_energy_landscape_plots(out_dir +"bi_only_refold",
                                energy_landscape_bidirectional_only_refolding)                                    
    # plot each unfolding/refolding pair, along with their velocities...
    kw_unfold = dict(style_data=dict(color='b',alpha=0.3))
    kw_fold = dict(style_data=dict(color='r',alpha=0.3))
    to_y = lambda x: x*1e12
    to_x = lambda x: x*1e9
    for i,(un,re) in enumerate(zip(unfolding_objs,refolding_objs)):
        # get the separation changes 
        min_v,max_v = min(un.Separation),max(un.Separation)
        fudge = (max_v-min_v)*0.1
        # using separation ('x value') as plotted y here, so use to_x
        ylim = to_x(np.array([min_v-fudge,max_v+fudge]))    
        fig = PlotUtilities.figure(figsize=(4,7))
        plt.subplot(2,1,1)
        FEC_Plot._fec_base_plot(x=un.Time,y=to_y(un.Force),**kw_unfold)
        FEC_Plot._fec_base_plot(x=re.Time,y=to_y(re.Force),**kw_fold)    
        PlotUtilities.lazyLabel("","Force","")        
        plt.subplot(2,1,2)
        FEC_Plot._fec_base_plot(x=un.Time,y=to_x(un.Separation),label="unfold",
                                **kw_unfold)
        FEC_Plot._fec_base_plot(x=re.Time,y=to_x(re.Separation),label="refold",
                                **kw_fold)  
        PlotUtilities.lazyLabel("","Separation","")        
        PlotUtilities.no_x_label()
        plt.ylim(ylim)        
        sep_unfold = to_x(un.ZFunc(un))
        plt.plot(un.Time,sep_unfold,color='k',linestyle='--')        
        plt.plot(re.Time,to_x(re.ZFunc(re)),label="Schedule",
                 color='k',linestyle=':')
        plt.ylim(ylim)
        PlotUtilities.lazyLabel("Time","Separation","")
        PlotUtilities.legend()
        name = out_dir + "unfolding_vs_time_{:d}.png".format(i)        
        PlotUtilities.savefig(fig,name)
    # make the plots we want                           
    for i,(ex,ex_unfold) in enumerate(zip(examples,unfolding_retracts)):
        name = out_dir + "{:d}".format(i)
        kw = dict(filter_fraction=1e-3)
        hairpin_plots(ex,out_path=name,**kw)
        # plot the unfolding force vs sep 
        fig = PlotUtilities.figure()
        FEC_Plot._fec_base_plot(x=ex_unfold.Separation,y=ex_unfold.Force)
        PlotUtilities.lazyLabel("Sep (nm)","Force (pN)","")
        PlotUtilities.savefig(fig,name + "_unfold.png")
        # plot the folding and refolding force vs time..
        fig = PlotUtilities.figure()
        FEC_Plot._fec_base_plot(x=unfolding_objs[i].Time,
                                y=to_y(unfolding_objs[i].Force))
        FEC_Plot._fec_base_plot(x=refolding_objs[i].Time,
                                y=to_y(refolding_objs[i].Force))                                
        PlotUtilities.lazyLabel("Time (s)","Force (pN)","")
        PlotUtilities.savefig(fig,name + "_bidirectional_fold.png")
    # XXX plot the data with the fit of the WLC aligned ontop 
        
if __name__ == "__main__":
    run()
