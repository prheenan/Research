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

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    abs_dir = "./"
    _, examples = FEC_Util.read_and_cache_pxp(abs_dir,force=False,
                                              cache_name="./cache.pkl")                    
    n_filter = 400                                          
    dicts = [dict(time_min=5.73,time_max=7.73),
             dict(time_min=5.99,time_max=7.17)]
    for i,ex in enumerate(examples):
        fec_pred,pred = Detector._predict_full(ex,threshold=1e-1,
                                               tau_fraction=0.005)
        fig = PlotUtilities.figure()
        Plotting.plot_prediction_info(fec_pred,pred)
        PlotUtilities.savefig(fig,"./pred{:d}.png".format(i))
        fig = PlotUtilities.figure()
        plt.subplot(2,1,1)
        FEC_Plot.force_versus_time(ex,NFilterPoints=n_filter)
        plt.subplot(2,1,2)
        FEC_Plot.z_sensor_vs_time(ex,NFilterPoints=n_filter)
        PlotUtilities.savefig(fig,"./out{:d}.png".format(i))
        # slice just the region we care about...
        ex = FEC_Util.slice_by_time(ex,**dicts[i])
        ex.Separation -= min(ex.Separation)
        ex.Force -= min(ex.Force)
        ex.Force *= -1
    # POST: examples are sliced
    # get the unfolding and refolding objects
    unfold,refold = [IWT_Util.get_unfold_and_refold_objects(e,number_of_pairs=1) 
                     for e in examples]
    # concatenate all of the refolding objects

    
if __name__ == "__main__":
    run()
