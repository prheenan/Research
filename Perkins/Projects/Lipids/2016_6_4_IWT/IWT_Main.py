# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D

sys.path.append("../../../../..")
# XXX move to utility?

from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util 
from GeneralUtil.python import PlotUtilities,CheckpointUtilities


def run():
    BaseDir = "/Volumes/group/4Patrick/Reports/" + \
           "2016_6_5_NearEquilibrium_Biotin_DOPC/IWT/"
    InBase = BaseDir + "In/"
    OutBase = BaseDir + "Out/"
    FullNames = [
    InBase +"2016-6-3-micah-1-part-per-million-biolevel-long-strept-coated.pxp",
    InBase +"2016-6-4-micah-1ppm-biolever-long-strept-saved-data.pxp",
    InBase +"2016-6-5-micah-1ppm-biolever-long-strept-saved-data.pxp"
    ]
    Limit = 100
    ForceReRead = False
    ForceRePlot = False
    IwtObjects,RetractList,Touchoff,LandscapeObj  = \
            CheckpointUtilities.getCheckpoint(OutBase + "IWT.pkl",
                                              IWT_Util.GetObjectsAndIWT,
                                              ForceReRead,InBase,FullNames,
                                              ForceReRead,
                                              Limit=Limit)
    ext,force  = CheckpointUtilities.\
        getCheckpoint(OutBase + "ExtAndForce.pkl",
                      IWT_Util.GetAllExtensionsAndForceAndPlot,
                      ForceRePlot,RetractList,
                      Touchoff,IwtObjects,OutBase)
    fig = PlotUtilities.figure(figsize=(8,12))
    IWT_Util.ForceExtensionHistograms(force,ext)
    PlotUtilities.savefig(fig,OutBase + "HeatMap.png")
    fig = PlotUtilities.figure(figsize=(8,12))
    IWT_Util.EnergyLandscapePlot(LandscapeObj)
    PlotUtilities.savefig(fig,OutBase + "IWT.png")

if __name__ == "__main__":
    run()
