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
from GeneralUtil.python import PlotUtilities
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from Research.Personal.EventDetection.Util import Analysis,Plotting 
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util



def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    data = FEC_Util.ReadInData("./Hold.pxp")
    split,info = Detector._predict_full(data[0],threshold=0.5)
    fig = PlotUtilities.figure(figsize=(4,8))
    Plotting.plot_prediction_info(split,info)
    PlotUtilities.savefig(fig,"./out.png")
    
    

if __name__ == "__main__":
    run()
