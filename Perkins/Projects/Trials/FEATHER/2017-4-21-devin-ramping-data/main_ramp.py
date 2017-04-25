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
from Research.Personal.EventDetection._2SplineEventDetector import Detector
from Research.Personal.EventDetection.Util import Analysis 
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util



def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    data = FEC_Util.ReadInData("./Graph3.pxp")
    split,info = Detector._predict_full(data[0],threshold=1e-1)
    plt.plot(split.retract.Force)
    for i in info.event_idx:
        plt.axvline(i)
    plt.show()

if __name__ == "__main__":
    run()
