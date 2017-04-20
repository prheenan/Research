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
from Research.Personal.EventDetection.Util import Analysis
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
    example_split,pred_info = Detector._predict_full(example,threshold=1e-4)
    retract = example_split.retract
    sep = example_split.retract.Separation
    sep -= min(sep)
    last = pred_info.event_idx[-1]
    last_x = sep[last]
    delta_x = 40e-9    
    force = example_split.retract.Force
    sep_bounds = [ [delta_x*1.2,last_x],
                   [0,delta_x]]
    models = CheckpointUtilities.getCheckpoint("./model{:d}.pkl".format(i),
                                               get_wlc_information,
                                               force_run,sep,force,sep_bounds,
                                               Ns=30)
    return models,retract                                                
        
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
        # get the FJC model...
        models,retract = CheckpointUtilities.getCheckpoint(
            "./model_all{:d}.pkl".format(i),get_basic_information,
            force_run,i,example)
        args.append([models,retract])
    models_all = [ [x0_x_y_tuple[0] for x0_x_y_tuple in list_v[0]] 
                   for list_v in args]
    first_l = np.concatenate([m[1] for m in models_all])
    last_l = np.concatenate([m[0] for m in models_all])
    print(first_l)
    style_common = dict(alpha=0.3)
    style_line = dict(linewidth=3,linestyle='--')
    fig = PlotUtilities.figure()                                                       
    plt.hist(first_l*1e9,label="first FJC",**style_common)
    plt.hist(last_l*1e9,label="second FJC",**style_common)
    plt.axvline(43,**style_line)
    plt.axvline(67,**style_line)
    PlotUtilities.lazyLabel("Contour Length (nm)","Count","")
    PlotUtilities.savefig(fig,"./out/o_hist{:d}.png".format(i))   
    for i,(models,retract) in enumerate(args):
        fig = PlotUtilities.figure()                                                   
        plt.plot(retract.Separation,retract.Force)
        for x0,model_x,model_y in models:
            plt.plot(model_x,model_y)
        plt.xlim([min(retract.Separation),max(retract.Separation)])
        PlotUtilities.lazyLabel("Separation (nm)","Force (pN)","")
        PlotUtilities.savefig(fig,"./out/out{:d}.png".format(i))
    
if __name__ == "__main__":
    run()
