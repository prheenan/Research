# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from GeneralUtil.python import GenUtilities
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = "//perknas2/group/4Patrick/CuratedData/Lipids/DOPC/"+\
            "NegativeControls/Representative_Gallery/"
    files = GenUtilities.getAllFiles(base,ext=".pxp")
    raw_data = IWT_Util.ReadInAllFiles(files,Limit=50)
    plt.plot(raw_data[0].Separation,RawData[0].Force)
    plt.show()

if __name__ == "__main__":
    run()
