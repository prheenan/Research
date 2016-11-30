# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from IgorUtil.PythonAdapter import PxpLoader

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    Base = "./"
    OutBase = Base + "out/"
    InFiles = [Base + "ForPatrick.pxp"]
    RawData = IWT_Util.\
              ReadInAllFiles(InFiles,Limit=50,
                             ValidFunc=PxpLoader.valid_fec_allow_endings)
    print(RawData)
    Example = RawData[0]
    plt.figure()
    plt.plot(Example.Separation,Example.Force)
    plt.show()
if __name__ == "__main__":
    run()
