# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,cProfile

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util import Analysis,Plotting
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
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
    # XXX debugging...
    time = retract.Time
    sep = example_split.retract.Separation
    sep -= min(sep)
    last = pred_info.event_idx[-1]
    last_x = sep[last]
    delta_x = 40e-9    
    force = example_split.retract.Force
    sep_bounds = [ [delta_x*1.2,last_x],
                   [0,delta_x]]
    """                            
    models = CheckpointUtilities.getCheckpoint("./model{:d}.pkl".format(i),
                                               get_wlc_information,
                                               force_run,sep,force,sep_bounds,
                                               Ns=30)
    """
    models = None
    return models,retract,pred_info                                               
        

def slice_retract(r,inf):
    n_points = int(0.02 * r.Separation.size)
    filtered = Analysis.filter_fec(r,n_points)
    idx,interp = Analysis.spline_interpolator_by_index(r.Force,n_points)
    where_force_above_zero = np.where(filtered.Force >= 0)[0]
    dt = r.Time[1] - r.Time[0]
    n_pairs =3
    t_per_pair= 2.2
    fudge = 0.5
    t_total = t_per_pair * n_pairs - fudge
    idx_to_move = int(np.ceil(t_total/dt))
    """
    plt.plot(r.Time[:last_le_sep_idx],r.Force[:last_le_sep_idx],
              color='k',alpha=0.3)
    plt.plot(r.Time[:last_le_sep_idx],interp(idx)[:last_le_sep_idx])
    plt.axvline(r.Time[last_idx_check[-1]])
    plt.show()
    """
    zero_offset = where_force_above_zero[0]
    slice_v = slice(zero_offset,zero_offset+idx_to_move)
    slice_obj = FEC_Util.MakeTimeSepForceFromSlice(r,slice_v)
    return slice_obj

def analyze_hairpin(abs_dir):
    _, examples = FEC_Util.read_and_cache_pxp(abs_dir,force=False)
    args = []
    force_run = False
    for i,example in enumerate(examples):
        # need to fix the dwell times; igor does not record it when using
        # the indenter
        example.set_dwell_time(example.Meta.DwellTime1)
        # get the FJC model...
        models,retract,pred_info = CheckpointUtilities.getCheckpoint(
            "./model_all{:d}.pkl".format(i),get_basic_information,
            force_run,i,example)
        args.append([models,retract,pred_info])
    # make a heat map of all the retracts...
    retracts = [a[1] for a in args]
    pred_info = [a[2] for a in args]
    just_ramping_portions = []
    for r,inf in zip(retracts,pred_info):
        just_ramp_tmp = slice_retract(r,inf)
        just_ramping_portions.append(just_ramp_tmp)
    # get split into unfolding and refolding regions...
    unfold,refold = [],[]
    for r in just_ramping_portions:
        unfold_tmp,refold_tmp = IWT_Util.\
            get_unfold_and_refold_objects_by_sep(r,number_of_pairs=3,
                                                 flip_forces=False,
                                                 fraction_for_vel=0.1)
        unfold.extend(unfold_tmp)
        refold.extend(refold_tmp)
        """
        for u,r in zip(unfold,refold):
            plt.subplot(2,1,1)
            plt.plot(u.Time,u.Force)
            plt.subplot(2,1,2)
            plt.plot(r.Time,r.Force)
            plt.show()
        """
    """

    Do a *really* fast-and-dirty energy landscape analysis
    """        
    keywords_iwt = dict(UnfoldingObjs=unfold,NumBins=75,RefoldingObjs=refold)
    landscape = CheckpointUtilities.\
        getCheckpoint("./landscape.pkl",
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
        
        y= ("Tilted (kT)(tilt: {:.2g}pN)".format(f))
        plt.ylim([-10,60])
        PlotUtilities.lazyLabel("Extension (nm)",y,"",frameon=True)
        ylim = np.array(plt.ylim())
        PlotUtilities.secondAxis(ax,label="kcal/mol",limits=ylim*0.529)        
        PlotUtilities.savefig(fig,"./out/landscape_{:.2g}.png".format(f*1e12))
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(just_ramping_portions,separation_max=100,
                           cmap='gist_earth')
    PlotUtilities.savefig(fig,"./out/heat.png")
    exit(1)
    x_plot = lambda x: x
    y_plot = lambda y: y*1e12
    n_filter_points = 500
    for i,(models,retract,pred_info) in enumerate(args):
        fig = PlotUtilities.figure(figsize=(8,12))        
        plt.subplot(2,1,1)
        FEC_Plot._fec_base_plot(retract.Separation * 1e9,
                                y_plot(retract.Force),
                                n_filter_points=n_filter_points)   
        PlotUtilities.lazyLabel("Separation (nm)","Force (pN)","")
        ax = plt.subplot(2,1,2)
        x_plot_tmp = x_plot(retract.Time)
        FEC_Plot._fec_base_plot(x_plot_tmp,
                                y_plot(retract.Force),
                                n_filter_points=n_filter_points)
        plt.xlim([min(x_plot_tmp),max(x_plot_tmp)])
        PlotUtilities.lazyLabel("Time (s)","Force (pN)","")
        retract_nm = 1e9 * retract.Separation
        limit_second = [min(retract_nm),max(retract_nm)]
        ax2 = PlotUtilities.secondAxis(ax,label="Separation (nm)",
                                       limits=limit_second,
                                       secondY=True,color='b')                                            
        ax2.plot(retract.Time,retract_nm,color='b')
        PlotUtilities.savefig(fig,"./out/out{:d}.png".format(i))    
    
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
    relative = [base_data_dir + "/3_ramps/",
                base_data_dir + "//"
    
        
    for abs_dir in [base_data_dir + r for r in relative]:
    analyze_hairpin(abs_dir)
    abs_dir = base_data_dir + \
              

    
if __name__ == "__main__":
    run()
