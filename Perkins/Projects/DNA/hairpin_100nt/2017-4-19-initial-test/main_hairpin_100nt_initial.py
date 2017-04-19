# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util import Analysis
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
from FitUtil.FreelyJointedChain.Python.Code import FJC

def get_single_example(abs_dir):
    _,force_extension_curves = FEC_Util.read_and_cache_pxp(abs_dir,
                                                           force=False)
    example = force_extension_curves[0]
    return example

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
    example = CheckpointUtilities.getCheckpoint("./single.pkl",get_single_example,
                                                False,abs_dir)
    example_split = Analysis.zero_and_split_force_extension_curve(example)
    sep = example_split.retract.Separation
    sep -= min(sep)
    force = example_split.retract.Force
    sep_bounds = [ [45e-9,71e-9],
                   [10e-9,40e-9]]
    models = []                   
    for min_sep,max_sep in sep_bounds:
        where_fit = np.where( (sep <= max_sep) & (sep >= min_sep))
        sep_fit = sep[where_fit]
        brute_dict = dict(ranges=[[max_sep/2,max_sep*1.5]],Ns=10)
        x0,model_x,model_y = FJC.fit_fjc_contour(sep_fit,
                                                 force[where_fit],
                                                 Lp=0.75e-9,
                                                 kbT=4.1e-21,
                                                 K0=800e-12,
                                                 brute_dict=brute_dict)
        models.append([x0,model_x,model_y]) 
    plt.plot(sep,force)
    for x0,model_x,model_y in models:
        print(x0)
        plt.plot(model_x,model_y)
    plt.show()
    
if __name__ == "__main__":
    run()
