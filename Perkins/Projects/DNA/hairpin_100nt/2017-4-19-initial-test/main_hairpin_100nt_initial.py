# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,cProfile,os

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util import Analysis,Plotting
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from GeneralUtil.python import PlotUtilities,CheckpointUtilities,GenUtilities
from FitUtil.FreelyJointedChain.Python.Code import FJC
from Research.Perkins.AnalysisUtil.EnergyLandscapes import \
   IWT_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass   
    
def get_wlc_information(sep,force,sep_bounds,**kwargs):
    models = []         
    for min_sep,max_sep in sep_bounds:
        where_fit = np.where( (sep <= max_sep) & (sep >= min_sep))
        sep_fit = sep[where_fit]
        min_range = np.mean([min_sep,max_sep])
        brute_dict = dict(ranges=[[0.9*max_sep,1.1*max_sep]],**kwargs)
        x0,model_x,model_y = FJC.fit_fjc_contour(sep_fit,
                                                 force[where_fit],
                                                 Lp=0.75e-9,
                                                 kbT=4.1e-21,
                                                 K0=800e-12,
                                                 brute_dict=brute_dict)
                    
        models.append([x0,model_x,model_y]) 
    return models
    
def get_basic_information(i,example,force_run=False):
    example_split = Analysis.zero_and_split_force_extension_curve(example)
    example_split,pred_info = Detector._predict_full(example,threshold=1e-1,
                                                     tau_fraction=0.02)
    retract = example_split.retract
    sep = example_split.retract.Separation
    sep -= min(sep)
    models = None
    return models,retract,pred_info                                               
        

def slice_retract(r,inf,n_pairs,slice_rel):
    n_points = int(0.02 * r.Separation.size)
    filtered = Analysis.filter_fec(r,n_points)
    idx,interp = Analysis.spline_interpolator_by_index(r.Force,n_points)
    where_force_above_zero = np.where(filtered.Force >= 0)[0]
    dt = r.Time[1] - r.Time[0]
    # XXX add in dt...
    idx_to_move = int(np.ceil((slice_rel.stop-slice_rel.start)/dt))
    zero_offset =where_force_above_zero[0]
    plt.plot(r.Time,r.Force,color='k',alpha=0.3)
    plt.plot(r.Time,interp(idx))
    plt.axvline(r.Time[zero_offset])
    plt.axvline(r.Time[min(r.Time.size-1,zero_offset+idx_to_move)])
    plt.show()
    slice_v = slice(zero_offset,zero_offset+idx_to_move)
    slice_obj = FEC_Util.MakeTimeSepForceFromSlice(r,slice_v)
    return slice_obj
   

def analyze_hairpin(abs_dir,force=False,**kw):
    id = "np={:d}_".format(kw["n_pairs"])
    cache_name = "./cache/cache{:s}.pkl".format(id)
    _, examples = FEC_Util.read_and_cache_pxp(abs_dir,force=force,
                                              cache_name=cache_name)
    args = []
    force_run = False
    dir_relative = os.path.dirname(abs_dir)
    rel = "./" + id
    GenUtilities.ensureDirExists("./out{:s}".format(id))                                                 
    for i,example in enumerate(examples[18:]):
        # get the FJC model...
        models,retract,pred_info = CheckpointUtilities.getCheckpoint(
            rel + "cache/model_all{:d}.pkl".format(i),get_basic_information,
            force_run,i,example)
        args.append([models,retract,pred_info])
    # make a heat map of all the retracts...
    retracts = [a[1] for a in args]
    pred_info = [a[2] for a in args]
    just_ramping_portions = []
    for r,inf in zip(retracts,pred_info):
        just_ramp_tmp = slice_retract(r,inf,**kw)
        just_ramping_portions.append(just_ramp_tmp)
    # get split into unfolding and refolding regions...
    unfold,refold = [],[]
    unfold_per_obj,refold_per_obj = [],[]
    for ramp in just_ramping_portions:
        try:
            unfold_tmp,refold_tmp = IWT_Util.\
                get_unfold_and_refold_objects_by_sep(ramp,
                                                     number_of_pairs=kw['n_pairs'],
                                                     flip_forces=False,
                                                     fraction_for_vel=0.1)
        except IndexError as e:
            print(e)
            continue
        unfold.extend(unfold_tmp)
        refold.extend(refold_tmp)
        unfold_per_obj.append(unfold_tmp)
        refold_per_obj.append(refold_tmp)
        """
        for u,r in zip(unfold_tmp,refold_tmp):
            plt.subplot(2,1,1)
            plt.plot(ramp.Time,ramp.Force,color='k',alpha=0.3)                
            plt.plot(u.Time,u.Force)
            plt.subplot(2,1,2)
            plt.plot(ramp.Time,ramp.Force,color='k',alpha=0.3)                            
            plt.plot(r.Time,r.Force)
            plt.show()
        """
 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(just_ramping_portions,separation_max=100,
                           cmap='gist_earth')
    PlotUtilities.savefig(fig,"./out{:s}/heat.png".format(id))
    x_plot = lambda x: x
    y_plot = lambda y: y*1e12
    n_filter_points = 500
    for i,retract in enumerate(just_ramping_portions):
        fig = PlotUtilities.figure(figsize=(8,12))        
        plt.subplot(4,1,1)
        FEC_Plot._fec_base_plot(retract.Separation * 1e9,
                                y_plot(retract.Force),
                                n_filter_points=n_filter_points)   
        PlotUtilities.lazyLabel("Separation (nm)","Force (pN)","")
        ax = plt.subplot(4,1,2)
        x_plot_tmp = x_plot(retract.Time)
        FEC_Plot._fec_base_plot(x_plot_tmp,
                                y_plot(retract.Force),
                                n_filter_points=n_filter_points)
        xlim =  [min(x_plot_tmp),max(x_plot_tmp)]                               
        plt.xlim(xlim)
        PlotUtilities.lazyLabel("","Force (pN)","")
        plt.subplot(4,1,3)        
        plt.plot(retract.Time,retract.Force*1e12,alpha=0.3)
        for u,r in zip(unfold_per_obj[i],refold_per_obj[i]):
            plt.plot(u.Time,u.Force*1e12,color='r',alpha=0.6)                
            plt.plot(r.Time,r.Force*1e12,color='b',alpha=0.6)   
        PlotUtilities.lazyLabel("","Force (pN)","")            
        plt.xlim(xlim)        
        plt.subplot(4,1,4)                         
        retract_nm = 1e9 * retract.Separation
        plt.plot(retract.Time,retract_nm,color='b')
        PlotUtilities.lazyLabel("Time (s)","Separation (nm)","")
        plt.xlim(xlim)                
        PlotUtilities.savefig(fig,"./out{:s}/{:d}.png".format(id,i))    
    """
    Do a *really* fast-and-dirty energy landscape analysis
    """        
    keywords_iwt = dict(UnfoldingObjs=unfold,NumBins=150,RefoldingObjs=refold)
    landscape = CheckpointUtilities.\
        getCheckpoint(rel + "landscape.pkl",
                      InverseWeierstrass.FreeEnergyAtZeroForce,
                      False,**keywords_iwt)
    L0_gc_poor_nm = 0.67 * (100 - 4 - 16*2)
    linker_lengths_nm = 3.2*2
    for f in [6e-12,8e-12,12e-12,16e-12]:
        tilted = IWT_Util.TiltedLandscape(landscape,f_one_half_N=f)   
        fig = PlotUtilities.figure()    
        ax = plt.subplot(2,1,1)    
        plt.plot(landscape.Extensions*1e9,landscape.EnergyLandscape/4.1e-21)
        ylim = np.array(plt.ylim())
        PlotUtilities.lazyLabel("","Free Landscape (kT)","")    
        PlotUtilities.secondAxis(ax,label="kcal/mol",limits=ylim*0.529)
        PlotUtilities.no_x_label()
        ax = plt.subplot(2,1,2)
        plt.plot(tilted.landscape_ext_nm,tilted.OffsetTilted_kT)
        plt.axvline(L0_gc_poor_nm,linestyle=':',linewidth=3,color='g',
                    label="L$_0$ (GC poor)")
        linker_kw = dict(color='r',linestyle='--',linewidth=3)
        plt.axvline(L0_gc_poor_nm+linker_lengths_nm,label="$\pm$ Linkers",
                    **linker_kw)
        plt.axvline(L0_gc_poor_nm-linker_lengths_nm,**linker_kw)
        
        y= ("Tilted (kT)(tilt: {:.2g}pN)".format(f*1e12))
        plt.ylim([-10,60])
        PlotUtilities.lazyLabel("Extension (nm)",y,"",frameon=True)
        ylim = np.array(plt.ylim())
        PlotUtilities.secondAxis(ax,label="kcal/mol",limits=ylim*0.529)        
        PlotUtilities.savefig(fig,"./out{:s}/landscape_{:.2g}.png".\
                              format(id,f*1e12))        
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base_data_dir = FEC_Util.default_data_root() + \
                    "4Patrick/CuratedData/DNA/hairpin-100nt-16gc/Positive"
    dict_10_ramps = dict(n_pairs =9,
                         slice_rel=slice(0.14,3.543))
    dict_3_ramps = dict(n_pairs =3,
                        slice_rel=slice(1.299,7.61))                
    dirs = [ [base_data_dir + "/3_ramps/",dict_3_ramps],
             [base_data_dir + "/10_ramps/",dict_10_ramps]]
    for abs_dir,kw in dirs:
        analyze_hairpin(abs_dir,**kw)
              

    
if __name__ == "__main__":
    run()
