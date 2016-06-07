# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


import copy
sys.path.append("../../../../..")
# XXX move to utility?
from IgorUtil.PythonAdapter import PxpLoader,ProcessSingleWave
from IgorUtil.PythonAdapter.WaveDataGroup import WaveDataGroup
from IgorUtil.PythonAdapter.DataObj import DataObj
from IgorUtil.PythonAdapter.TimeSepForceObj import TimeSepForceObj,Bunch

from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python.IgorUtil import SavitskyFilter
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
from GeneralUtil.python.LibUtil.CypherUtil import ConvertSepForceToZsnsrDeflV

def ReadInData(FullName,Limit=3):
    """
    Reads in the PXP waves as TimeSepForce Object (must *only* wave waves with
    Cypher-like wave endins in FullName
    
    Args:
        FullName: Path to the pxp file to load in
        Limit: maximum number to read in. 
    """
    MData = PxpLoader.LoadPxp(FullName)
    # convert the waves into TimeSepForce objects
    Objs = [TimeSepForceObj(WaveDataGroup(v)) for _,v in MData.items()]
    return Objs[:Limit]

def ToPullingObjects(TimeSepForceObjects):
    """
    Converts TimeSepForceObjects to InverseWeierstrass objects

    Args:
        TimeSepForceObjects: list of TimeSepForceObjects to transform
    """
    Objs = [InverseWeierstrass.\
            FEC_Pulling_Object(Time=o.Time,
                               Extension=o.Separation,
                               Force=o.Force,
                               SpringConstant=o.SpringConstant,
                               Velocity=o.Velocity)
            for o in TimeSepForceObjects]
    return Objs



def MakeTimeSepForceFromSlice(Obj,Slice):
    """
    Given a TimeSepForceObject and a slice, gets a new object using just the slices 
    given

    Args:
        Obj:
        Slice:
    """
    ToRet = TimeSepForceObj()
    GetSlice = lambda x: x[Slice].copy()
    ToRet.LowResData = DataObj(GetSlice(Obj.Time),
                               GetSlice(Obj.Separation),
                               GetSlice(Obj.Force),
                               copy.deepcopy(Obj.Meta))
    return ToRet

def BreakUpIntoApproachAndRetract(mObjs):
    """
    Takes in a list of TimeSepForceObj, returns a list of approach and retract 
    objects, where index [i] in both lists refers to original curve i.

    Args:
        mObjs: Lsit of TimeSepForce Objects
    Returns:
        Tuple of Two lists: Approach,Retract, which are the TimeSepForce
        Objects of the approach and retract, respectively
    """
    Approach = []
    Retract = []
    for o in mObjs:
        ForceArray = o.Force
        SepArray = o.Separation
        # note: force is 'upside down' by default, so high force (near surface
        # is actually high)
        MinForceIdx = np.argmax(ForceArray)
        # get the different slices
        SliceAppr = slice(0,MinForceIdx)
        SliceRetr = slice(MinForceIdx,None)
        # Make a new object with the given force and separation
        # at approach and retract
        Approach.append(MakeTimeSepForceFromSlice(o,SliceAppr))
        Retract.append(MakeTimeSepForceFromSlice(o,SliceRetr))
    return Approach,Retract

def GetAroundTouchoff(Objects,fraction=0.05,FilterPoints=20,
                      MetersAfterTouchoff=None):
    """
    Gets the data 'MetersAfterTouchoff' meters after (in ZSnsr)the surface touchoff, 
    based on taking the median of
    'fraction' points far away from the surface

    Args:
        Objects: list of TimeSepForce Objects
        fraction: Amount to average to determine the zero point for the force. 
        FilterPoints: how many points to filter to find the zero 
        MetersFromTouchoff: gets this many meters away from the surface. If
        None, just returns all the data
    Returns:
        List of TimeSepForce objects, each offset in force and separation to 
        zero at the 
    """
    ToRet = []
    for o in Objects:
        ForceArray = o.Force
        SepArray = o.Separation
        ForceFilter = SavitskyFilter(o.Force,nSmooth=FilterPoints)
        # Flip the sign of the force
        ForceSign = -1 * ForceFilter 
        MinForceIdx = np.argmin(ForceSign)
        N = ForceSign.size
        NumMed = int(N*fraction)
        MedRetr = np.median(ForceSign[-NumMed:])
        ZeroForce = ForceSign - MedRetr
        # Get the first time the Retract forces is above zero
        FilteredRetract = SavitskyFilter(ZeroForce)
        ZeroIdx = np.where(FilteredRetract >= 0)[0][0] + MinForceIdx
        ZSnsr,_ = o.ZsnsrAndDeflV
        # Get just that part of the Retract
        StartRetractX = ZSnsr[ZeroIdx]
        if (MetersAfterTouchoff is not None):
            EndRetractX = StartRetractX + MetersAfterTouchoff
            Index = np.arange(0,N)
            # XXX build in approach/retract
            StopIdxArr = np.where( (ZSnsr > EndRetractX) &
                                   (Index > ZeroIdx))[0][0]
        else:
            # just get eveything
            StopIdxArr = None
        NewSlice = slice(ZeroIdx,StopIdxArr)
        MyObj = MakeTimeSepForceFromSlice(o,NewSlice)
        ToRet.append(MyObj)
    return ToRet

def run():
    Base = "./Data/"
    FullName = Base + "2016-6-5-micah-1ppm-biolever-long-strept-saved-data.pxp"
    mObjs = pCheckUtil.getCheckpoint(Base + "cache.pkl",ReadInData,False,
                                     FullName)
    ApproachList,RetractList = BreakUpIntoApproachAndRetract(mObjs)
    PastZeroExt = 8e-9
    FilterPoints = 30
    # get just after the touchoff
    Touchoff = GetAroundTouchoff(RetractList,FilterPoints=30,
                                 MetersAfterTouchoff=PastZeroExt)
    fig = plt.figure()
    for Retract,Touch in zip(RetractList,Touchoff):
        ForZ = Retract
        RetractZ,_ = Retract.ZsnsrAndDeflV
        RetractZ -= np.min(RetractZ)
        plt.subplot(2,1,1)
        plt.plot(RetractZ,Retract.Force,alpha=0.3)
        plt.subplot(2,1,2)
        TouchZ,_ = Touch.ZsnsrAndDeflV
        # sign correct and offset the force
        Force = Touch.Force * -1
        Force -= Force[0]
        # Zero out the ZSnsr
        TouchZ -= min(TouchZ)
        # Zero out the Sep
        Sep = Touch.Separation
        Sep -= Sep[0]
        plt.plot(Touch.Separation,Force,alpha=0.3)
        fig.savefig("out.png")

if __name__ == "__main__":
    run()
