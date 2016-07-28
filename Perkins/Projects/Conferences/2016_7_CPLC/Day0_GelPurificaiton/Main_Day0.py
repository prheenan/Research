# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../../")
from GeneralUtil.python import PlotUtilities as pPlotUtil
from FitUtil.FitUtils.Python import FitUtil as FitUtil
from Research.Perkins.AnalysisUtil.Gels import ImageJUtil



def run():
    """
    Creates a printout giving the mean and stdev of the 
    concentraiton measurements
    """
    # all concentrations are (control,1 minute digestion)
    # these concentrations as in ng/uL
    ConcentrationsNanodrop = [33,20.5]
    # next, we read in the gel
    LaneObj = ImageJUtil.ReadFileToLaneObj("ControlAndDigestionLane.xls")
    # Get the two concentrations; control is first
    ConcentrationsGels = LaneObj.Lanes
    # okay, get everything on a scale relative to the control
    Relative = lambda x: np.array(x)/x[0]
    # we *only* want the digestion lane, or element 1
    RelativeDigestions = [Relative(ConcentrationsNanodrop)[1],
                          Relative(ConcentrationsGels)[1]]
    MeanRel = np.mean(RelativeDigestions)
    StdevRel = np.std(RelativeDigestions)
    print("The mean post-1 minute digestion content relative to the undigested"+
          " control is {:.2f}+/-{:.2f}".format(MeanRel,StdevRel))
    print("This implies {:.2f}+/- {:.2f} of the DNA wasn't able to circularize"\
          .format(1-MeanRel,StdevRel))
    
if __name__ == "__main__":
    run()
