# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import GeneralUtil.python.GenUtilities as pGenUtil
from collections import OrderedDict
sys.path.append("../../../../")

class LaneObject:
    def __init__(self,*lanes):
        self.Lanes = np.array(list(lanes))
        self.TotalIntensity = sum(self.Lanes)
    def Normalized(self,LaneIndex):
        assert LaneIndex < len(self.Lanes) , "Asked for a lane we dont have"
        return self.Lanes[LaneIndex]/self.TotalIntensity
    def __str__(self):
        return "\n".join("Lane{:03d}={:.2f}".format(i,self.Normalized(i))
                         for i in range(self.Lanes.size))

def GetImageJData(DataDirBase,ext=".xls"):
    """
    Given a base data directory, finds all files with ext in each subdirectory

    Args:
        DataDirBase: base data directory. Each subdirectory has files with 
        extension 'ext'
       
        ext: file extension
    Returns:
        ordered dictionary of <subdir:fullpaths>
    """
    Voltages = OrderedDict()
    for f in sorted(os.listdir(DataDirBase)):
        PossibleSubDir = DataDirBase + f +"/"
        if (os.path.isdir(PossibleSubDir)):
            Files = pGenUtil.getAllFiles(PossibleSubDir,".xls")
            Voltages[f] =Files
    return Voltages

def ReadFileToLaneObj(File):
    return LaneObject(GetImageJMeasurements(File))

def GetImageJMeasurements(File):
    """
    Returns the in-order values of the intensity column in the ImageJ xls file

    Args:
        File: to read from
    Returns:
        intensity column
    """
    return np.loadtxt(File,skiprows=1,usecols=(1,))
