# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from GeneralUtil.python import GenUtilities,CheckpointUtilities
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
    base_debug = "../../_1ReadDataToCache/debug_no_event/"
    out_dir = "./out/"
    files = GenUtilities.getAllFiles(base_debug,ext=".pkl")
    data = [CheckpointUtilities.getCheckpoint(f,None,False) for f in files]
    for fec in data:
        split_fec = Analysis.zero_and_split_force_extension_curve(fec)
        files = GenUtilities.ensureDirExists(out_dir)
        _,info = Detector._predict_full(fec,threshold=1e-1)
        # XXX fix threshhold
        fig = PlotUtilities.figure(figsize=(8,12))    
        Plotting.plot_prediction_info(split_fec,info)
        PlotUtilities.savefig(fig,"{:s}{:s}.png".format(out_dir,fec.Meta.Name))



if __name__ == "__main__":
    run()
