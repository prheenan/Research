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
    where_force_above_zero = np.where(filtered.Force >= 0)[0]
    # the slice of the retract is between the zero point and the last time 
    # we are (effectively) at the surface
    last_event_idx = inf.event_idx[-1]
    min_sep = np.percentile(filtered.Separation[:last_event_idx],50)    
    where_sep_le_min_sep = np.where(filtered.Separation <= min_sep)[0]
    slice_v = slice(where_force_above_zero[0],where_sep_le_min_sep[-1])
    slice_obj = FEC_Util.MakeTimeSepForceFromSlice(r,slice_v)
    return slice_obj

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base_data_dir = FEC_Util.default_data_root()
    abs_dir = base_data_dir + \
              "4Patrick/CuratedData/DNA/hairpin-100nt-16gc/Positive"
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
        if (i == 5):
            break
    # make a heat map of all the retracts...
    retracts = [a[1] for a in args]
    pred_info = [a[2] for a in args]
    just_ramping_portions = []
    for r,inf in zip(retracts,pred_info):
        just_ramp_tmp = slice_retract(r,inf)
        just_ramping_portions.append(just_ramp_tmp)
        plt.plot(just_ramp_tmp.Force)
        for i in range(0,3):
            plt.axvline((just_ramp_tmp.size/3) * i)
        plt.show()
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(just_ramping_portions,separation_max=100,
                           cmap='gist_earth')
    PlotUtilities.savefig(fig,"./out/heat.png")
    exit(1)
    """
    Do a *really* fast-and-dirty energy landscape analysis
    """
    """
    """
    models_all = [ [x0_x_y_tuple[0] for x0_x_y_tuple in list_v[0]] 
                   for list_v in args]
    first_l = np.concatenate([m[1] for m in models_all])
    last_l = np.concatenate([m[0] for m in models_all])
    style_common = dict(alpha=0.3)
    style_line = dict(linewidth=3,linestyle='--')
    fig = PlotUtilities.figure()                                                       
    plt.hist(first_l*1e9,label="first FJC",**style_common)
    plt.hist(last_l*1e9,label="second FJC",**style_common)
    plt.axvline(43,**style_line)
    plt.axvline(67,**style_line)
    PlotUtilities.lazyLabel("Contour Length (nm)","Count","")
    PlotUtilities.savefig(fig,"./out/o_hist{:d}.png".format(i))   
    """
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
    """
    
if __name__ == "__main__":
    run()
