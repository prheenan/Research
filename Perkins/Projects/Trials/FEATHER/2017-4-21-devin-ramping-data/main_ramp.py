# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,scipy
sys.path.append("../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from Research.Personal.EventDetection.Util import Analysis,Plotting,InputOutput
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass,WeierstrassUtil
from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import UnivariateSpline


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    file_path = "./2017-6-5-devin-email-Hold.pxpImage1056Concat.csv"
    data = InputOutput.read_and_cache_file(file_path=file_path,
                                           cache_directory="./",
                                           has_events=True,force=False)
    split,info = Detector._predict_full(data,threshold=0.5)
    fig = PlotUtilities.figure(figsize=(4,8))
    Plotting.plot_prediction_info(split,info)
    PlotUtilities.savefig(fig,"./out.png")
    
    

if __name__ == "__main__":
    run()
