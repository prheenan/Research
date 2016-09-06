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
    Gets the volume to dilute a sample with a given volume and concentration
    
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
        if (len(Array) > 1):
            return Array[idx]
        else:
            return Array[0]
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


def SeriallyDilute(Stock,DesiredConcentrations,DesiredVolumes):
    """
    Given a stock and desired concentraitons and desired volumes at each 
    concentration, returns the list of stocks, volumes, dilutions, and final
    stocks

    Args:
         Stock: concentration, same units as elements of DesiredConcentrations
         DesiredConcentrations: array or scalar, same units as stock. These
         are the concentrations we want

         DesiredVolumes: scalar or array volumes, in units of au/<Stock>,
         we want for each dilution. Note that *actual* volume will be a bit 
         more, since we need something to serially dilute with. E.g., if 
         DesiredVolumes was 1L, we might need 10mL extra for 'downstream'
         Dilutions
    Returns
        Tuple of arrays, the elements are grouped from high to low 
        concentrations for each of:<What stocks we used, What volumes of stocks,
        what volume we diluted with, what was the resulting stock>.
        Note that the first and last elements are just Stock and DesiredVolumes
    """
    NumConc = len(DesiredConcentrations)
    MassNeededBelow = 0
    VolumeStock = []
    VolumeDilute = []
    Stocks = []
    ResultingStock = []
    # work backwards with the last one first, determine what volume
    # and concentraton is needed
    for i in range(NumConc-1,-1,-1):
        VolumeNeeded = GetFromArrayOrScalar(DesiredVolumes,i)
        ConcNeeded = GetFromArrayOrScalar(DesiredConcentrations,i)
        # what mass is needed 'below' us?
        MassNeeded = ConcNeeded*VolumeNeeded
        MassNeededBelow += MassNeeded
        # we use  the stock 'above' what we need here
        TmpStock = Stock if (i==0) \
                   else GetFromArrayOrScalar(DesiredConcentrations,i-1)
        # determine how much stock we need
        VolStock = MassNeededBelow/TmpStock
        # how can we dilute it to get the concentration we need?
        # since we have already accounting for the underlying mass,
        # this implicitly includes the 'extra' volume for dilutions
        VolDilute = GetVolumeToDilute(TmpStock,VolStock,
                                      ConcNeeded)
        VolumeStock.append(VolStock)
        VolumeDilute.append(VolDilute)
        Stocks.append(TmpStock)
        ResultingStock.append(ConcNeeded)
    # reverse the arrays so we go big to small (natural order for dilution)
    RetSanitize = lambda x: x[::-1]
    RetArrs = Stocks,VolumeStock,VolumeDilute,ResultingStock
    return [RetSanitize(a) for a in RetArrs]

def PrintSerialDilution(Stocks,VolumeStock,VolumeDilute,FinalStocks,
                        VolString="uL",ConcString="ng/uL"):
    """
    Given the results of SeriallyDilute, prints off the relevant information
    to 

    Args:
        Stocks,VolumeStock,VolumeDilute,FinalStocks: output of SeriallyDilute
        VolString,ConcStrung:  units for the volume and concentration
    """
    for stock,VolStock,VolDilute,DilutedStock in \
        zip(Stocks,VolumeStock,VolumeDilute,FinalStocks):
        VolumeTotal = VolStock + VolDilute
        StockStr = "{:5.3g}{:s} of {:5.3g}{:s} with {:5.3g}{:s} Buffer".\
                   format(VolStock,VolString,stock,ConcString,VolDilute,
                          VolString)
        ResultStr = "{:5.3g}{:s} of {:5.3g}{:s}".\
                    format(VolumeTotal,VolString,DilutedStock,ConcString)
        print("{:s} gives {:s}".format(StockStr,ResultStr))

def PrintSerialSteps(Stock,Volumes,Desired,
                     ConcString="ng/uL",VolString="uL"):
    """
    Given a stock concentration, desired volumes and concentrations, prints
    out the steps needed to serially dilute

    Args:
        see PrintSerialDilution
    """
    Stocks,VolumeStock,VolumeDilute,FinalStocks = \
                        SeriallyDilute(Stock,Desired,Volumes)
    PrintSerialDilution(Stocks,VolumeStock,VolumeDilute,
                        FinalStocks,ConcString=ConcString,
                        VolString=VolString)

def PrintSolutionSteps(Stats,Volume,vol_units="uL"):
    """
    Prints the steps to seriall dilute things
    Args:
        Stats: List of Tuples; each element is <Name,Concentration Unit,
        Stock Concentraiton, Desired concentration, mass present in solution
        already> 
    """
    # get the stocks, desired concntrations, and already-present concentraitons
    Stocks = [s[2] for s in Stats]
    Desired = [s[3] for s in Stats]
    Already = [s[4] for s in Stats]
    Volumes = GetVolumesNeededByConcentration(Stocks,Desired,Volume,
                                              AlreadyHaveMass=Already)
    BufferVolume = Volume - sum(Volumes)
    # check that our buffer is reasonable non-negative. if it is very close
    # to zero (less than 1% error), let it slide.
    assert (BufferVolume > -Volume/100) , \
        "Warning: cant make this solution. Need a negative volume of buffer. "+\
        "Use more concentrated stocks"
    print("In a total solution of {:.1f}{:s}...".format(Volume,vol_units))
    for (name,conc_units,conc_stock,desired_conc,_),vol_stock in\
        zip(Stats,Volumes):
        print("\t{:.2f}{:s} of {:.2f}{:s} {:s} for {:.2f}{:s} in solution".\
              format(vol_stock,vol_units,conc_stock,conc_units,name,
                     desired_conc,conc_units))
    print("\tRemainder is ({:.1f}{:s}) of buffer".format(BufferVolume,
                                                         vol_units))
  

    

