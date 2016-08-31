# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np

class DilutionObj:
    def __init__(self,StockConc,StockVol,DesiredConc,AddVol,Name=""):
        self.StockConc=StockConc
        self.StockVol=StockVol
        self.DesiredConc=DesiredConc
        self.AddVol=AddVol
        self.Name = Name

class SolutionObj:
    def __init__(self,ArrayOfDilutionObjects):
        self.arr = ArrayOfDilutionObjects

def GetVolumeToDilute(Concentration,Volume,DesiredConcentration):
    """
    Gets the volume to dilute a sampel with a given volume and concentration
    
    Args:
        Concentration: mass per volume of the system
        Volume: volume of the system.  Units of mass/Concentration
        DesiredConcentration: desired concentration after the diltion
    Returns:
        Amount of additional volume to add to the system; total volume of
        Volume + (<return of this funciton>) gives the desired concentration
    """
    # figure out how much stuff we have
    ng = Concentration*Volume
    # what total volume do we need?
    volumeNeeded = ng/DesiredConcentration
    # how much do we need to add?
    volumeToAdd = volumeNeeded-Volume
    return volumeToAdd

def PrintDilutions(StockConcs,StockVols,ConcDesired,UnitVol=None,
                   UnitConc=None,**kwargs):
    """
    Convenience wrapper: Print all the dilutions, given desired concentrations
    stocks, etc.

    Args:
        StockConcs,StockVols,ConcDesired:  see GetDilutionObj
        UnitConc,UnitVol: see PrintVolumeDilutions
        **kwargs: ppassed directly to GetDilutionObj
    Returns:
        all the dilution objects
    """
    if (UnitVol is None):
        UnitVol = ["uL" for _ in StockConcs]
    if (UnitConc is None):
        UnitConc = ["ng/uL" for _ in StockConcs]    
    dilutionObj = GetDilutionObjects(StockConcs,StockVols,ConcDesired,**kwargs)
    _PrintVolumeDilutions(dilutionObj,UnitVol=UnitVol,UnitConc=UnitConc)
    return dilutionObj
    
def GetDilutionObjects(StockConcs,StockVols,ConcDesired,**kwargs):
    """
    Args: 
        StockConcs,StockVols,ConcDesired: see GetDilutionObj
        **kwargs: ppassed directly to GetDilutionObj
    Returns:
        List of dilution objects
    """
    dilutions = GetVolumeToDilute(StockConcs,StockVols,ConcDesired)
    dilutionObj = [GetDilutionObj(StockConcs,StockVols,ConcDesired,d,i,
                                  **kwargs)
                   for i,d in enumerate(dilutions)]
    return dilutionObj

def GetFromArrayOrScalar(Array,idx):
    """
    Tries to get Array[idx]; otherwise just returns Array
    
    Args:
        Array: array-like
        idx: number 
    Returns:
        relevant element of the array
    """
    try:
        return Array[idx]
    except (TypeError,IndexError) as e:
        return Array
    
def GetDilutionObj(StockConcs,StockVols,ConcDesired,VolToAdd,Idx,
                   StockName=""):
    """
    Returns a Dilution object at a given index, given all the informaiton
    
    Args:
        StockConcs: array or scalar-like of stock concentrations
        StockVols: array or scalar-like of stock volumes
        ConcDesired: array or scalar-like of desired concentrations
        VolToAdd: array or scalar-like of volumes to add
        Idx: index within all the arrays we want
        StockName: optional names of the stock
    Returns:
        DilutionObj to use
    """
    DesiredConc = GetFromArrayOrScalar(ConcDesired,Idx)
    StockVol = GetFromArrayOrScalar(StockVols,Idx)
    StockConc = GetFromArrayOrScalar(StockConcs,Idx)
    AddVol =  GetFromArrayOrScalar(VolToAdd,Idx)
    StockName = GetFromArrayOrScalar(StockName,Idx)
    return DilutionObj(StockConc=StockConc,StockVol=StockVol,
                       DesiredConc=DesiredConc,
                       AddVol=AddVol,Name=StockName)

def _DilutionString(dilutionObj,UnitVol,UnitConc):
    """
    Args:
        dilutionObj: list of dilution objects
        UnitVol: string, unit of volume
        UnitConc: string, unit of Concentration
    Returns:
        String representation of dilution objects
    """
    toRet = ""
    for i,d in enumerate(dilutionObj):
        stockConcStr = "({:4.1f}{:s}@{:4.1f}{:s})".\
                       format(float(d.StockVol),UnitVol[i],float(d.StockConc),
                              UnitConc[i])
        volAddStr = "({:4.1f}{:s})".format(d.AddVol,UnitVol[i])
        TotalVol = float(d.AddVol) + float(d.StockVol)
        toRet += ("Stock {: <7} (#{:03d}) {: <20}" +
                  "add {: <10} -> {:3.1f}{:7s} in {:3.1f}{:7s}\n").\
            format(d.Name,i,stockConcStr,volAddStr,d.DesiredConc,
                   UnitConc[i],TotalVol,UnitVol[i])
    return toRet
    
def _PrintVolumeDilutions(dilutionObj,**kwargs):
    """
    Gets and prints all the dilution objects

    Args:
        dilutionObj: list of dilution objects
        **kwargs: see DilutionString
    """
    print(_DilutionString(dilutionObj,**kwargs))
    
def GetVolumesNeededByConcentration(StockConcs,ConcsDesired,TotalVolume,
                                    AlreadyHaveMass=None):
    """
    Given desired and stock concentrations and a final volume, gives the 
    volumes needed of each stock

    Args:
        StockConcs: index i refers to some species, same units as 
        ConcsDesired[i]
     
        ConcsDesired: what we want in the volume

        AlreadyHaveMass: if present, the mass already present in the buffer
        we will use. Element [i] should have the same units as StockConcs[i]
    Returns:
        Array of volumes needed going from StockConcs to ConcsDesired in 
        TotalVolume (note that TotalVolume-sum(<Return of this function>) is
        taken up by some unspecified buffer)
    """
    if (AlreadyHaveMass is None):
        AlreadyHaveMass  = np.zeros_like(StockConcs)
    StockArr = np.array(StockConcs)
    TotalVolumeNeeded = np.array(ConcsDesired)*TotalVolume/StockArr
    EffectiveVolumeAlreadyPresent = np.array(AlreadyHaveMass)/StockArr
    return TotalVolumeNeeded - EffectiveVolumeAlreadyPresent

    

