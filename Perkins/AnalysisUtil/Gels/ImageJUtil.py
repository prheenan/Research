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


class LaneObject(object):
    def __init__(self,*lanes):
        self.Lanes = np.array(list(lanes))
        self.TotalIntensity = np.sum(self.Lanes)
    def Normalized(self,LaneIndex):
        assert LaneIndex < len(self.Lanes) , "Asked for a lane we dont have"
        return self.Lanes[LaneIndex]/self.TotalIntensity
    def __str__(self):
        return "\n".join("Lane{:03d}={:.2f}".format(i,self.Normalized(i))
                         for i in range(self.Lanes.size))

class OverhangLane(LaneObject):
    """
    Class to keep track of the bands in an overhang lane
    """
    def __init__(self,Linear,Circular=0,Concat=0):
        super(OverhangLane,self).__init__(Linear,Circular,Concat)
        self.LinearBand=Linear
        self.CircularBand = Circular
        self.Concatemers = Concat
    def _Norm(self,x):
        return x/self.TotalIntensity
    @property
    def LinearRelative(self):
        return self._Norm(self.LinearBand)
    @property
    def CircularRelative(self):
        return self._Norm(self.CircularBand)
    @property
    def ConcatemerRelative(self):
        return self._Norm(self.Concatemers)
    def __str__(self):
        return "Lin:{:3.2f},Circ:{:3.2f},Concat:{:3.2f}".\
            format(self.LinearRelative,
                   self.CircularRelative,
                   self.ConcatemerRelative)
    def __repr__(self):
        return str(self)
    

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
    """
    Given a file, get it as a lane object.

    Args:
        File: see ReadFileToOverhangObj
    Returns:
        LaneObject
    """
    return LaneObject(*list(GetImageJMeasurements(File)))

def ReadFileToOverhangObj(File):
    """
    get the file's data, as an OverhangLane object

    Args:
       File: full path to the file to read. Must be formatted with the second
       column being the intensitieis, and the rows being <concat if any,
       circular if any, linear>
    Returns:
       OverhangLane object
    """
    # reverse so it goes linear,circular,conat
    Measurements = list(GetImageJMeasurements(File)[::-1])
    return OverhangLane(*Measurements)
    
def GetImageJMeasurements(File):
    """
    Returns the in-order values of the intensity column in the ImageJ xls file

    Args:
        File: to read from
    Returns:
        intensity column
    """
    return np.loadtxt(File,skiprows=1,usecols=(1,))
