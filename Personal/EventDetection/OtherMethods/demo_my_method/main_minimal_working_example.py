# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from Research.Personal.EventDetection.OtherMethods import method_helper
from Research.Personal.EventDetection.Util import Analysis,Plotting,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector

from GeneralUtil.python import PlotUtilities

    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    fec = method_helper.get_example()
    thresh = 2e-2
    split_fec = Analysis.zero_and_split_force_extension_curve(fec)   
    _,info = Detector._predict_full(fec,threshold=thresh)
    # XXX fix threshhold
    fig = PlotUtilities.figure(figsize=(8,12))    
    Plotting.plot_prediction_info(split_fec,info)
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
