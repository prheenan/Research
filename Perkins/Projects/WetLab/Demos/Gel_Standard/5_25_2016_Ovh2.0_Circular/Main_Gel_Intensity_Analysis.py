# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import ntpath

class GelObject:
    def __init__(self,BaseDir,**kwargs):
        self.BaseDir = BaseDir
        for k,v in kwargs.items():
            setattr(self,k,v)
        self._Lanes = []
        self._Names = []
    def AddLanes(self,ToSetData,ToSetNames):
        self._Lanes.append(ToSetData)
        self._Names.append(ToSetNames)

def GetLogFile(Dir):
    """
    Given a directory, return the log file per the conventions
    
    Args:
        Dir : the directory of interest
    Returns:
        String to the log file
    """
    return Dir + "log.txt"

def GetLanesFolder(Dir,ImageJName="ImageJ"):
    return Dir + ImageJName + "/" 

def GetLanesFiles(Dir,ext=".xls"):
    return [GetLanesFolder(Dir) + f
            for f in os.listdir(GetLanesFolder(Dir)) if f.endswith(ext)]

def HasImageJFolder(Dir):
    return os.path.isdir(GetLanesFolder(Dir)) and \
        os.path.isfile(GetLogFile(Dir))

def GetDataObject(Dir):
    # function to remove bad characers etc
    SanitizeName = lambda s: s.split("(")[0].replace(" ","").replace("%","")
    # get the log (meta) file
    with open(GetLogFile(Dir)) as f:
        fields =[line.split(":") for line in f.read().split("\n")]
        fields = [ (SanitizeName(f[0]),f[1]) for f in fields if len(f[0]) > 0]
        metadict = dict(fields)
        toRet = GelObject(Dir,**metadict)
    # get the lanes files
    Lanes = GetLanesFiles(Dir)
    Delim = "\t"
    ExtSize = 4
    for LaneFile in Lanes:
        # get the data (first column is just numbering)
        LaneData = np.loadtxt(LaneFile,skiprows=1,delimiter=Delim)[:,1]
        NumLanes = LaneData.size
        # get the file name
        FileName = ntpath.basename(LaneFile)
        # Remove the extension
        FileName = FileName[:-ExtSize]
        # split the file name to get the labels
        Labels = [o for o in FileName.split("_")]
        # just get the last N, where N is the number of lanes
        Labels = Labels[-NumLanes:]
        toRet.AddLanes(LaneData,Labels)
    return toRet


    

def Plot(DataObjects):
    Voltages = [o.Voltage for o in DataObjects]
    # get the normalized circular yields for linear and circular
    CircularName = "circular"
    LinearName = "linear"
    CircularAbsArr = []
    LinearAbsArr = []
    TotalAbsArr = []
    VoltagesArr = []
    for k,o in enumerate(DataObjects):
        Names = [[n.lower() for n in listV] for listV in o._Names]
        LaneIdxWithCircular = [i for i,lane in enumerate(o._Lanes)
                               if CircularName.lower() in Names[i]]
        CircularAbs = [band
                       for j in LaneIdxWithCircular
                       for k,band in enumerate(o._Lanes[j])
                       if CircularName.lower()== Names[j][k]]
        LinearAbs = [band
                     for j in LaneIdxWithCircular
                     for k,band in enumerate(o._Lanes[j])
                     if LinearName.lower()== Names[j][k]]
        TotalAbs = [sum(o._Lanes[j]) for k in LaneIdxWithCircular]
        CircularAbsArr.extend(CircularAbs)
        LinearAbsArr.extend(LinearAbs)
        TotalAbsArr.extend(TotalAbs)
        VoltagesArr.extend([float(o.Voltage) for _ in CircularAbs])
    # convert everything to arrays
    CircularAbs = np.array(CircularAbsArr)
    LinearAbs = np.array(LinearAbsArr)
    TotalAbs = np.array(TotalAbsArr)
    Voltages = np.array(VoltagesArr)
    MinV = min(Voltages)
    MaxV = max(Voltages)
    Range = MaxV-MinV
    # fudge is for plotting: how much to extend plot in either x direction.
    Fudge = Range/10
    xlim = [MinV-Fudge,MaxV+Fudge]
    plt.plot(Voltages,CircularAbs/(LinearAbs+CircularAbs),'ro')
    plt.xlabel("Voltages")
    plt.ylabel("Band Intensity of Circular DNA, Relative to Total Yield")
    plt.xlim(xlim)
    plt.show()
    

def AnalyzeInDirectory(Dir):
    SubDirs = [Dir + o + "/"
               for o in os.listdir(Dir) if os.path.isdir(Dir + o) ]
    # get only those subdirectories matching our function
    Func = HasImageJFolder
    SubDirsMatching = [o for o in SubDirs if Func(o)]
    # Get the load files and associated data files
    DataObjects= [] 
    for d in SubDirsMatching:
        DataObjects.append(GetDataObject(d))
    #Plot them all
    Plot(DataObjects)
    
              
def run():
    """
    Runs gel analysis on the given directory
    """
    Dir = "/Users/patrickheenan/Documents/education/boulder_files/"+\
          "rotations_year_1/3_perkins/data_gels_nanodrop/gel/5_2016_summer/"
    AnalyzeInDirectory(Dir)

if __name__ == "__main__":
    run()
