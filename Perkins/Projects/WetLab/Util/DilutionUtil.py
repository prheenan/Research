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
    try:
        ng = Concentration*Volume
    except TypeError:
        ng = np.array(Concentration)*np.array(Volume)
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
        UnitVol = ["uL" for _ in ConcDesired]
    if (UnitConc is None):
        UnitConc = ["ng/uL" for _ in ConcDesired]    
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
    n = len(dilutionObj)
    for i,d in enumerate(dilutionObj):
        stockConcStr = "({:4.1f}{:s}@{:4.1f}{:s})".\
                       format(float(d.StockVol),UnitVol[i],float(d.StockConc),
                              UnitConc[i])
        volAddStr = "({:4.1f}{:s})".format(d.AddVol,UnitVol[i])
        TotalVol = float(d.AddVol) + float(d.StockVol)
        toRet += ("{: <4} (#{:03d}) {: <20}" +
                  "add {: <8} -> {:3.1f}{:6s} in {:3.1f}{:7s}").\
            format(d.Name,i,stockConcStr,volAddStr,d.DesiredConc,
                   UnitConc[i],TotalVol,UnitVol[i])
        if (i != n-1):
            toRet += "\n"
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
        we will use. Element [i] should have the same 'mass' units as 
        StockConcs[i]
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


def SeriallyDilute(Stock,DesiredConcentrations,DesiredVolumes,
                   dilution_concentration=0):
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

         dilution_concentration: the concentration of whatever already in the 
         stock.  (i.e. if we aren't using something with a concentration of 
         zero. For example, if diluting 100mM NaCl with 10mM dilution,
         Stock would be 100, DilutionConcentration would be 10
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
        # determine what total volume of the final solution we need
        # (we know the mass, and the concentration is specified) 
        V0 = MassNeededBelow/ConcNeeded
        TmpStock = Stock if (i==0) \
                   else GetFromArrayOrScalar(DesiredConcentrations,i-1)
        conc_diff = dilution_concentration - TmpStock
        # We are solving the following system:
        # c_stock * V_s + c_dilute * V_dilute = MassNeededBelow
        # V_s + V_dilute                      = V0
        VolStock = (dilution_concentration*V0-MassNeededBelow)/conc_diff
        VolDilute = (MassNeededBelow-TmpStock*V0 )/conc_diff
        # we use  the stock 'above' what we need here
        VolumeStock.append(VolStock)
        VolumeDilute.append(VolDilute)
        Stocks.append(TmpStock)
        ResultingStock.append(ConcNeeded)
    # reverse the arrays so we go big to small (natural order for dilution)
    RetSanitize = lambda x: x[::-1]
    RetArrs = Stocks,VolumeStock,VolumeDilute,ResultingStock
    return [RetSanitize(a) for a in RetArrs]

def PrintSerialDilution(Stocks,VolumeStock,VolumeDilute,FinalStocks,
                        VolString="uL",ConcString="ng/uL",
                        BufferString="Buffer"):
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
        StockStr = "{:5.3g}{:s} of {:5.3g}{:s} with {:5.3g}{:s} {:s}".\
                   format(VolStock,VolString,stock,ConcString,VolDilute,
                          VolString,BufferString)
        ResultStr = "{:5.3g}{:s} of {:5.3g}{:s}".\
                    format(VolumeTotal,VolString,DilutedStock,ConcString)
        print("{:s} gives {:s}".format(StockStr,ResultStr))

def StockVolumeNeededForSerialDilution(Stock,Volumes,Desired):
    """
    Gets the  total volume needed of the 'stock'

    Args:
         see PrintSerialSteps
    """
    _,VolumeStock,_,_ = SeriallyDilute(Stock,Desired,Volumes)
    return VolumeStock[0]

        
def PrintSerialSteps(Stock,Volumes,Desired,
                     ConcString="ng/uL",VolString="uL",BufferString="Buffer",
                     **kwargs):
    """
    Given a stock concentration, desired volumes and concentrations, prints
    out the steps needed to serially dilute

    Args:
        see PrintSerialDilution
    """
    Stocks,VolumeStock,VolumeDilute,FinalStocks = \
                        SeriallyDilute(Stock,Desired,Volumes,**kwargs)
    PrintSerialDilution(Stocks,VolumeStock,VolumeDilute,
                        FinalStocks,ConcString=ConcString,
                        VolString=VolString,BufferString=BufferString)

def PrintSolutionSteps(Stats,Volume,vol_units="uL",BufferName="buffer",
                       PostVolume=0):
    """
    Prints the steps to seriall dilute things
    Args:
        Stats: List of Tuples; each element is <Name,Concentration Unit,
        Stock Concentraiton, Desired concentration, mass present in solution
        already> 
        PostVolume: if true, this is the volume to add after some step (e.g.
        thawing). We store the solution at a higher concentration
    """
    # get the stocks, desired concntrations, and already-present concentraitons
    Stocks = [s[2] for s in Stats]
    Desired = [s[3] for s in Stats]
    Already = [s[4] for s in Stats]
    Volumes = GetVolumesNeededByConcentration(Stocks,Desired,Volume,
                                              AlreadyHaveMass=Already)
    BufferVolume = Volume - sum(Volumes) - PostVolume
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
    print("\tRemainder is ({:.1f}{:s}) of {:s}".format(BufferVolume,
                                                       vol_units,BufferName))
    if (PostVolume > 1e-12):
        print("\tTo use, add ({:.1f}{:s}) of {:s}".format(PostVolume,
                                                          vol_units,BufferName))

  

    

