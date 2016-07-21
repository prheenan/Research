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

class DNAWlcPoints:
    class BoundingIdx:
        def __init__(self,start,end):
            self.start = start
            self.end = end
    def __init__(self,Wlc1Start,Wlc1End,Wlc2Start,Wlc2End,TimeSepForceObj):
        self.FirstWLC = DNAWlcPoints.BoundingIdx(Wlc1Start,Wlc1End)
        self.SecondWLC = DNAWlcPoints.BoundingIdx(Wlc2Start,Wlc2End)
        self.TimeSepForceObject = TimeSepForceObj

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


def CalculateOffset(Obj,Slice):
    """
    Calculates the force Offset for the given object and slice
    """
    return np.median(Obj.Force[Slice])

def ZeroForceAndSeparation(Obj,IsApproach,FractionForOffset=0.1):
    """
    Given an object and whether it is the approach or retract, zeros it
    out 

    Args:
        Obj: TimeSepForce Object
        IsApproach: True if this represents the approach portion of the 
        Curve.
    Returns:
        Tuple of <Zeroed Separation, Zeroed Force>
    """
    Separation = Obj.Separation
    Time = Obj.Time
    SeparationZeroed = Separation - min(Separation)
    N = SeparationZeroed.size
    NOffsetPoints = int(np.ceil(N))
    # approach's zero is at the beginning (far from the surface)
    # retract is at the end (ibid)
    if (IsApproach):
        SliceV = slice(0,NOffsetPoints,1)
    else:
        SliceV = slice(-NOffsetPoints,None,1)
    # get the actual offset assocaiated with this object
    Offset = CalculateOffset(Obj,SliceV)
    return SeparationZeroed.copy(),(Obj.Force - Offset).copy()

def PreProcessApproachAndRetract(Approach,Retract,
                                 NFilterPoints=100,
                                 ZeroForceFraction=0.2,
                                 ZeroSep=True,FlipY=True):
    """
    Given and already-split approach and retract curve, pre=processes it.
    This *modifies the original array* (ie: in-place)

    Args:
        Approach,Retract: output of GetApproachRetract 
        NFilterPoints: number of points for finding the surface
        ZeroForceFraction: if not None, fraction of points near the retract end
        to filter to
        
        ZeroSep: if true, zeros the separation to its minima
        FlipY: if true, multiplies Y (force) by -1 before plotting
    Returns:
        tuple of <Appr,Retr>, both pre-processed TimeSepFOrce objects for
        the appropriate reason
    """
    if (ZeroForceFraction is not None):
        # then we need to offset the force
        # XXX assume offset is the same for both
        _,ZeroForceRetr = GetSurfaceIndexAndForce(Retract,
                                                  Fraction=ZeroForceFraction,
                                                  FilterPoints=NFilterPoints,
                                                  ZeroAtStart=False)
        _,ZeroForceAppr = GetSurfaceIndexAndForce(Approach,
                                                  Fraction=ZeroForceFraction,
                                                  FilterPoints=NFilterPoints,
                                                  ZeroAtStart=True)
        # add, because the sign diffreent presummably hasnt been fixed
        # (See below)
        Approach.Force += ZeroForceRetr
        # Do the same for retract
        Retract.Force += ZeroForceRetr
    if (ZeroSep):
        MinSep = min(np.min(Approach.Separation),
                     np.min(Retract.Separation))
        Approach.Separation -= MinSep
        Retract.Separation -= MinSep
    if (FlipY):
        Approach.Force *= -1
        Retract.Force *= -1
    return Approach,Retract
    
def PreProcessFEC(TimeSepForceObject,**kwargs):
    """
    Returns the pre-processed (zeroed, flipped, etc) approach and retract
    
    Args:
        TimeSepForceObject: the object we are dealing with. copied, not changed
        **kwargs: passed directly to PreProcessApproachAndRetract
    Returns: 
        tuple of pre-processed approach and retract, see 
        PreProcessApproachAndRetract
    """
    Appr,Retr = GetApproachRetract(TimeSepForceObject)
    # now pre-process and overwrite them
    Appr,Retr = PreProcessApproachAndRetract(Appr,Retr)
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
    Get the approach and retraction curves of a TimeSepForceObject. Does *not*
    include the dwell portion
    
    Args:
        o: the TimeSepForce Object, assumed 'raw' (ie: invols peak at top)
    Returns:
        TUple of <Appr,Retract>, which are both TimeSepForce object of the
        Approach and Retract regions
    """
    ForceArray = o.Force
    TimeEndOfApproach = o.TriggerTime
    TimeStartOfRetract = TimeEndOfApproach + o.SurfaceDwellTime
    # figure out where the indices we want are
    Time = o.Time
    IdxEndOfApproach = np.argmin(np.abs(Time-TimeEndOfApproach))
    IdxStartOfRetract = np.argmin(np.abs(Time-TimeStartOfRetract))
    # note: force is 'upside down' by default, so high force (near surface
    # is actually high) is what we are looking for
    # get the different slices
    SliceAppr = slice(0,IdxEndOfApproach)
    SliceRetr = slice(IdxStartOfRetract,None)
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
                            ZeroAtStart=True,FlipSign=True):
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

        FlipSign: if true (default), assumes the data is 'raw', so that
        Dwell happens at positive force. Set to false if already fixed
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
    if (FlipSign):
        ForceSign = -1 * ForceFilter
    else:
        ForceSign = ForceFilter
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
                        MetersAfterTouchoff=None,Correct=False,**kwargs):
    """
    Args:
        fraction: Amount to average to determine the zero point for the force. 
        FilterPoints: how many points to filter to find the zero, from the 
        *start* of the array forward
        MetersFromTouchoff: gets this many meters away from the surface. If
        None, just returns all the data
        Correct: if true, corrects the data by flipping it and zeroing it out. 

        **kwargs: passed on to GetSurfaceIndexAndForce
    """
    ZeroIdx,MedRetr =  GetSurfaceIndexAndForce(o,fraction,FilterPoints,
                                               ZeroAtStart=False,**kwargs)
    if (MetersAfterTouchoff is not None):
        XToUse  = o.Separation
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
    if (Correct):
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

def GetFECPullingRegionAlreadyFlipped(Retract,NFilterFraction=0.02):
    """
    Adapter to GetFECPullingRegion, gets the retract after, assuming things
    have already been flipped
    """
    # now, get just the 'post touchoff' region. We *dont* want to flip
    # the sign when doing this
    SurfaceNPoints = int(np.ceil(Retract.Force.size *NFilterFraction))
    Retract = GetFECPullingRegion(Retract,FlipSign=False,
                                  FilterPoints=SurfaceNPoints)
    return Retract

def FilteredGradient(Retract,NFilterFraction):
    """
    Get the filtered gradient of an object

    The force is filtered, the gradient is taken, the gradient is filtered

    Args:
       Retract: TimeSepForce Objectr
       NFilterFraction: see GetGradientOutliersAndNormalsAfterAdhesion

    Returns:
       Array of same size as Retract.force which is the fitlered gradient.
    """
    N = Retract.Force.size
    NFilterPoints = int(np.ceil(NFilterFraction*N))
    RetractZeroSeparation = Retract.Separation
    RetractZeroForce = Retract.Force
    FilteredForce = GetFilteredForce(Retract,NFilterPoints)
    FilteredForceGradient = np.gradient(FilteredForce.Force)
    return FilteredForceGradient

def GetGradientOutliersAndNormalsAfterAdhesion(Retract,NFilterFraction,
                                               NoTriggerDistance=150e-9):
    """
    Returns where we are past the notrigger distance, where we are outlying
    and where we are normal, according to a filtered first derivative

    Args:
        Retract: TimeSepForce Object to use, assumed zeroed
        NFilterFraction: [0,1], how much to use for Savitsky Filtering
        NoTriggerDistance: how long (in meters) not to trigger
        after reaching the surface
    Returns:
        Tuple of <Indices after adhesion, Indices outlying q75+iqr,Indices
         within q75>
    """
    RetractZeroSeparation = Retract.Separation
    FilteredForceGradient = FilteredGradient(Retract,NFilterFraction)
    NoAdhesionMask = RetractZeroSeparation > NoTriggerDistance
    OnlyPositive = np.abs(FilteredForceGradient[NoAdhesionMask])
    q75, q25 = np.percentile(OnlyPositive, [75 ,25])
    iqr = q75-q25
    # q75 + 1.5 * iqr is a good metric for outliers
    IsOutlier = lambda x: x > q75 + 1.5 * iqr
    # being less than q75 is a good sign we are normal-ish
    IsNormal = lambda x:  (x <= q75)
    # where does the first WLC start?
    Outliers = np.where(IsOutlier(FilteredForceGradient) & \
                        NoAdhesionMask)
    # determine where the normal, non adhesion points are
    NormalPoints = np.where( IsNormal(FilteredForceGradient) &
                             NoAdhesionMask )
    return NoAdhesionMask,Outliers[0],NormalPoints[0]

def GetWLCPoints(NoAdhesionMask,Outliers,Normal,Object):
    """
    Gets the indices associated with (near) the start and end of the WLC
    
    Args:
       see GetWlcIdxObject
    Returns:
        tuple of <start idx 1, end idx 1, start idx 2, end idx 2> where
        the numbers refer to the WLC region
    """
    FirstBackToNormalPost = lambda Point: Normal[np.where( (Normal > Point))][0]
    # first WLC starts with the first outlier
    FirstWlcStart = Outliers[0]
    # first WLC ends after we are back to 'normal'
    EndOfFirstWLC = FirstBackToNormalPost(FirstWlcStart)
    # now we have a (very rough) guess for the contour length; the next
    # WLC should be at least 70% later by XXX (literature). Say 50% to be safe
    NextWLCSizeMin = 1.5
    Separation = Object.Separation
    SeparationFirstWLC = Separation[EndOfFirstWLC]
    MinSeparationSecond = SeparationFirstWLC * NextWLCSizeMin
    # second WLC starts when we are an outlier again
    StartOfSecondWLC = Outliers[np.where( (Outliers > EndOfFirstWLC) &
                                          (Separation[Outliers]
                                           > MinSeparationSecond))][0]
    # ends when we ae back to normal 
    EndOfSecondWLC = FirstBackToNormalPost(StartOfSecondWLC)
    return FirstWlcStart,EndOfFirstWLC,StartOfSecondWLC,EndOfSecondWLC

def GetWlcIdxObject(NoAdhesionMask,Outliers,Normal,Retract):
    """
    Returns the DNAWlcPoints associated with the object

    Args:
        Retract: the rertract object used
        others: See output of GetGradientOutliersAndNormalsAfterAdhesion
    Returns:
        DNAWlcPoints object associated
    """
    Points = GetWLCPoints(NoAdhesionMask,Outliers,Normal,Retract)
    return DNAWlcPoints(*Points,TimeSepForceObj=Retract)
    
def GetRegionForWLCFit(RetractOriginal,NFilterFraction=0.1,**kwargs):
    """
    Given a (pre-processed, so properly 'flipped' and zeroed) WLC, gets the 
    approximate region for WLC fitting (ie: roughly the first curve up to 
    but not including the overstretching transition

    Args:
        RetractOriginal: pre-processed retract
        **kwargs: passed to GetFECPullingRegionAlreadyFlipped
    Returns:
        TimeSepForce Object of the portion of the curve to fit 
    """
    Retract = GetFECPullingRegionAlreadyFlipped(RetractOriginal,
                                                NFilterFraction=NFilterFraction,
                                                **kwargs)
    RetractZeroSeparation = Retract.Separation
    FilteredForceGradient = FilteredGradient(Retract,
                                             NFilterFraction=NFilterFraction)
    NoAdhesionMask,Outliers,NormalPoints = \
        GetGradientOutliersAndNormalsAfterAdhesion(Retract,
                                                   NFilterFraction)
    Idx = GetWlcIdxObject(NoAdhesionMask,Outliers,NormalPoints,Retract)
    FirstOutlier = Idx.FirstWLC.start
    EndOfFirstWLC = Idx.FirstWLC.end
    N = RetractOriginal.Force.size
    # determine the maximum point between the two outliers; this is likely
    # the linear (ie: stretching) point of the WLC
    MiddleOfFirstWLC = FirstOutlier + \
                    np.argmax(FilteredForceGradient[FirstOutlier:EndOfFirstWLC])
    IndexForWLCFit = int(np.mean([MiddleOfFirstWLC,EndOfFirstWLC]))
    MetersAfterTouchoff = RetractZeroSeparation[IndexForWLCFit]
    # *dont* sign correct anything.
    NearSurface = GetFECPullingRegion(Retract,
                                      MetersAfterTouchoff=MetersAfterTouchoff,
                                      Correct=False,FlipSign=False)
    return NearSurface

    

