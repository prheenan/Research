# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

import copy
from IgorUtil.PythonAdapter.DataObj import DataObj
from IgorUtil.PythonAdapter.TimeSepForceObj import TimeSepForceObj,Bunch
from IgorUtil.PythonAdapter import PxpLoader,ProcessSingleWave
from IgorUtil.PythonAdapter.WaveDataGroup import WaveDataGroup
from IgorUtil.PythonAdapter.TimeSepForceObj import TimeSepForceObj,Bunch
from GeneralUtil.python.IgorUtil import SavitskyFilter


def ReadInData(FullName,Limit=None):
    """
    Reads in the PXP waves as TimeSepForce Object (must *only* wave waves with
    Cypher-like wave endins in FullName
    
    Args:
        FullName: Path to the pxp file to load in
        Limit: maximum number to read in. If none, defaults to all (!)
    """
    MData = PxpLoader.LoadPxp(FullName)
    # convert the waves into TimeSepForce objects
    Objs = [TimeSepForceObj(WaveDataGroup(v)) for _,v in MData.items()]
    if (Limit is not None):
        return Objs[:Limit]
    else:
        return Objs

def MakeTimeSepForceFromSlice(Obj,Slice):
    """
    Given a TimeSepForceObject and a slice, gets a new object using just the 
    slices given
    Args:
        Obj:
        Slice:
    """
    ToRet = TimeSepForceObj()
    # note: we make a copy, to avoid any reference funny business
    GetSlice = lambda x: x[Slice].copy()
    ToRet.LowResData = DataObj(GetSlice(Obj.Time),
                               GetSlice(Obj.Separation),
                               GetSlice(Obj.Force),
                               copy.deepcopy(Obj.Meta))
    return ToRet


def UnitConvert(TimeSepForceObj,
                ConvertX=lambda x : x,
                ConvertY=lambda y : y,
                GetX = lambda x : x.Separation,
                GetY = lambda x : x.Force):
    """
    Converts the 'X' and 'Y' using the specified units and properties 
    of the object passed in 

    Args:
        TimeSepForceObj : see ApproachRetractCurve
        ConvertX: method to convert the X values into whatever units we want
        ConvertY: metohod to convery the Y values into whatever units we want
        GetX: gets the x values (assumed separation for plotting XXX TODO)
        GetY: gets the y values (assumed force for plotting XXX TODO)
    Returns: 
        deep *copy* of original object in the specified units
    """
    ObjCopy = copy.deepcopy(TimeSepForceObj)
    ObjCopy.Force = ConvertY(GetY(ObjCopy))
    ObjCopy.Separation = ConvertX(GetX(ObjCopy))
    return ObjCopy


def PreProcessFEC(TimeSepForceObject,NFilterPoints=100,
                  ZeroForceFraction=0.2,
                  ZeroSep=True,FlipY=True):
    """
    Returns the pre-processed (zeroed, flipped, etc) approach and retract
    
    Args:
        TimeSepForceObject: the object we are dealing with. copied, not changed
        NFilterPoints: number of points for finding the surface
        ZeroForceFraction: if not None, fraction of points near the retract end
        to filter to
        
        ZeroSep: if true, zeros the separation to its minima
        FlipY: if true, multiplies Y (force) by -1 before plotting
    Returns:
        tuple of <Appr,Retr>, both pre-processed TimeSepFOrce objects for
        the appropriate reason
    """
    Appr,Retr = GetApproachRetract(TimeSepForceObject)
    if (ZeroForceFraction is not None):
        # then we need to offset the force
        # XXX assume offset is the same for both
        _,ZeroForceRetr = GetSurfaceIndexAndForce(Retr,
                                                  Fraction=ZeroForceFraction,
                                                  FilterPoints=NFilterPoints,
                                                  ZeroAtStart=False)
        _,ZeroForceAppr = GetSurfaceIndexAndForce(Appr,
                                                  Fraction=ZeroForceFraction,
                                                  FilterPoints=NFilterPoints,
                                                  ZeroAtStart=True)
        Appr.Force += ZeroForceAppr
        # Do the same for retract
        Retr.Force += ZeroForceAppr
    if (ZeroSep):
        MinSep = np.min(TimeSepForceObject.Separation)
        Appr.Separation -= MinSep
        Retr.Separation -= MinSep
    if (FlipY):
        Appr.Force *= -1
        Retr.Force *= -1
    return Appr,Retr

def SplitAndProcess(TimeSepForceObj,ConversionOpts=dict(),
                    NFilterPoints=100,**kwargs):
    """
    Args:
        TimeSepForceObj: see PreProcessFEC
        ConversionOpts: passed to UnitConvert
        NFilterPoints: see PreProcessFEC
        **kwargs: passed to PreProcessFEC
    """
    # convert the x and y to sensible units
    ObjCopy = UnitConvert(TimeSepForceObj,**ConversionOpts)
    # pre-process (to, for example, flip the axes and zero everything out
    Appr,Retr = PreProcessFEC(ObjCopy,NFilterPoints=NFilterPoints,**kwargs)
    return Appr,Retr


def GetApproachRetract(o):
    """
    Get the approach and retraction curves of a TimeSepForceObject
    
    Args:
        o: the TimeSepForce Object, assumed 'raw' (ie: invols peak at top)
    Returns:
        TUple of <Appr,Retract>, which are both TimeSepForce object of the
        Approach and Retract regions
    """
    ForceArray = o.Force
    # note: force is 'upside down' by default, so high force (near surface
    # is actually high) is what we are looking for
    MinForceIdx = np.argmax(ForceArray)
    # get the different slices
    SliceAppr = slice(0,MinForceIdx)
    SliceRetr = slice(MinForceIdx,None)
    # Make a new object with the given force and separation
    # at approach and retract
    Appr = MakeTimeSepForceFromSlice(o,SliceAppr)
    Retr = MakeTimeSepForceFromSlice(o,SliceRetr)
    return Appr,Retr


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
        Appr,Retr = GetApproachRetract(o)
        Approach.append(Appr)
        Retract.append(Retr)
    return Approach,Retract

def GetFilteredForce(Obj,NFilterPoints):
    """
    Given a TimeSepForce object, return a (filtered) copy of it

    Args:
        Obj: the TimeSepForce object we care about
        NFilterPoitns: fed to savitsky golay filter
    """
    ToRet = copy.deepcopy(Obj)
    ToRet.Force = SavitskyFilter(Obj.Force,nSmooth=NFilterPoints)
    return ToRet

def GetSurfaceIndexAndForce(TimeSepForceObj,Fraction,FilterPoints,
                            ZeroAtStart=True):
    """
    Given a retraction curve, a fraction of end-points to take the median of,
    and a filtering for the entire curve, determines the surface location

    Args:
        TimeSepForceObj: single TimeSepForce Object, assumes 'raw', so that
        invols region is a negative force.
        Fraction: see GetAroundTouchoff
        FilterPoints: see GetAroundTouchoff
        ZeroAtStart: if true, uses the first 'fraction' points; otherwise 
        uses the last 'fraction' points
    Returns: 
        Tuple of (Integer surface index,Zero Force)
    """
    o = TimeSepForceObj
    ForceArray = o.Force
    SepArray = o.Separation
    if (FilterPoints > 1):
        ForceFilter = SavitskyFilter(o.Force,nSmooth=FilterPoints)
    else:
        ForceFilter = o.Force
    # Flip the sign of the force
    ForceSign = -1 * ForceFilter 
    N = ForceSign.size
    NumMed = int(N*Fraction)
    if (ZeroAtStart):
        SliceMed = slice(0,NumMed,1)
    else:
        SliceMed = slice(-NumMed,None,1)
    MedRetr = np.median(ForceSign[SliceMed])
    ZeroForce = ForceSign - MedRetr
    # Get the first time the Retract forces is above zero
    FilteredRetract = SavitskyFilter(ZeroForce)
    ZeroIdx = np.where(FilteredRetract >= 0)[0][0]
    return ZeroIdx,MedRetr

def GetFECPullingRegion(o,fraction=0.05,FilterPoints=20,
                        MetersAfterTouchoff=None):
    """
    fraction: Amount to average to determine the zero point for the force. 
    FilterPoints: how many points to filter to find the zero, from the 
    *start* of the array forward

    MetersFromTouchoff: gets this many meters away from the surface. If
    None, just returns all the data
    """
    ZeroIdx,MedRetr =  GetSurfaceIndexAndForce(o,fraction,FilterPoints,
                                               ZeroAtStart=False)
    if (MetersAfterTouchoff is not None):
        XToUse  = o.Zsnsr
        N = XToUse.size
        # Get just that part of the Retract
        StartRetractX = XToUse[ZeroIdx]
        EndRetractX = StartRetractX + MetersAfterTouchoff
        Index = np.arange(0,N)
        # XXX build in approach/retract
        StopIdxArr = np.where( (XToUse > EndRetractX) &
                               (Index > ZeroIdx))[0][0]
    else:
        # just get eveything
        StopIdxArr = None
    NewSlice = slice(ZeroIdx,StopIdxArr)
    MyObj = MakeTimeSepForceFromSlice(o,NewSlice)
    # sign correct and offset the force
    MyObj.Force = MyObj.Force * -1
    MyObj.Force -= MedRetr
    MyObj.Separation -= MyObj.Separation[0]
    MyObj.Zsnsr -= MyObj.Zsnsr[0]
    return MyObj
def GetAroundTouchoff(Objects,**kwargs):
    """
    Gets the data 'MetersAfterTouchoff' meters after (in ZSnsr)the surface 
    touchoff,  based on taking the median of
    'fraction' points far away from the surface

    XXX: Generalize to approach and retract

    Args:
        Objects: list of TimeSepForce Objects
        **kwargs: passed directly to GetPullingRegion
    Returns:
        List of TimeSepForce objects, each offset in force and separation to 
        zero at the perceived start
    """
    ToRet = []
    for o in Objects:
        ToRet.append(GetFECPullingRegion(o,**kwargs))
    return ToRet

