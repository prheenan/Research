# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")

from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import PlotUtilities

from scipy.signal import sawtooth
from scipy.optimize import fminbound

def RobTimeSepForceToIWT(o,ZFunc):
    k = o.Meta.__dict__["K"]
    velocity = o.Meta.__dict__["RetractVelocity"]
    Obj = InverseWeierstrass.FEC_Pulling_Object(Time=o.Time,
                                                Extension=o.Separation,
                                                Force=o.Force,
                                                SpringConstant=k,
                                                Velocity=velocity,
                                                ZFunc=ZFunc)
    Obj.SetWork(Obj.CalculateForceCummulativeWork())
    return Obj

def DistanceToRoot(DeltaA,Beta,ForwardWork,ReverseWork):
    """
    Gives the distance to the root in equation 18 (see NumericallyGetDeltaA)

    Args:
        Beta,DeltaA: See ibid
        FowardWork,ReverseWork: list of the works as defined in ibid, same
        Units as DeltaA
    """
    nf = len(ForwardWork)
    nr = len(ReverseWork)
    # get the forward and reverse 'factor': difference should be zero
    Forward = np.mean(1/(nr + nf * np.exp(Beta * (ForwardWork-DeltaA))))
    Reverse = np.mean(1/(nf + nr * np.exp(Beta * (ReverseWork+DeltaA))))
    # we really only case about the abolute value of the expression, since
    # we want the two sides to be equal...
    return np.abs(Forward-Reverse)

def NumericallyGetDeltaA(Forward,Reverse,disp=3,**kwargs):
    """
    Numerically solves for DeltaA, as in equation 18 of 

    Hummer, G. & Szabo, A. 
    Free energy profiles from single-molecule pulling experiments. 
    PNAS 107, 21441-21446 (2010).

    Note that we use a root finder to find the difference in units of kT,
    then convert back (avoiding large floating point problems associated with
    1e-21)

    Args:
        Forward: List of forward paths
        Reverse: List of reverse paths
    Returns:
        Free energy different, in joules
    """
    # XXX should fix/ not hard code
    beta = 1/(4.1e-21)
    # get the work in terms of beta, should make it easier to converge
    Fwd = [f.Work*beta for f in Forward]
    Rev = [f.Work*beta for f in Reverse]
    MaxWorks = [np.max(np.abs(Fwd)),
                np.max(np.abs(Rev))]
    MinWorks = [np.min(Fwd),
                np.min(Rev)]
    Max = max(MaxWorks)
    Min = min(MinWorks)
    # only look between +/- the max. Note that range is guarenteed positive
    Range = Max-Min
    FMinArgs = dict(x1=-Range,x2=Range,full_output=True,disp=disp,**kwargs)
    # note we set beta to one, since it is easier to solve in units of kT
    ToMin = lambda A: DistanceToRoot(A,Beta=1,ForwardWork=Fwd,ReverseWork=Rev)
    xopt,fval,ierr,nfunc = fminbound(ToMin,**FMinArgs)
    return xopt/beta
    
def run():
    """
    """
    Base = "./"
    OutBase = Base + "out/"
    InFiles = [Base + "RNA_Ramp_20nmps.pxp"]
    RawData = IWT_Util.ReadInAllFiles(InFiles,Limit=50)
    FilterPoints = 100
    FilteredData = [FEC_Util.GetFilteredForce(r,FilterPoints)
                    for r in RawData]
    Example = FilteredData[0]
    z0_Separation = 60e-9
    zf_Separation = 73e-9
    z0_ZSnsr = 57e-9
    zf_ZSnsr = 70e-9
    # choose the bounds
    z0 = z0_Separation
    zf = zf_Separation
    cycles = 16
    size = Example.ZSnsr.size
    cyclesize = int(size/cycles)
    retracts,reverse = [],[]
    for i in range(cycles):
        # XXX only get first half
        half = int((i+0.5)*cyclesize)
        fwd_slice = slice(i*cyclesize,half)
        rev_slice = slice(half,(i+1)*cyclesize)
        fwd = FEC_Util.MakeTimeSepForceFromSlice(Example,fwd_slice)
        rev = FEC_Util.MakeTimeSepForceFromSlice(Example,rev_slice)
        retracts.append(fwd)
        reverse.append(rev)
    fwd_cycle_size= fwd.Force.size
    rev_cycle_size = rev.Force.size
    up = np.linspace(z0,zf,fwd_cycle_size)
    down = np.linspace(zf,z0,rev_cycle_size)
    updown = np.concatenate((up,down))
    cat_cyc = np.concatenate([updown for _ in range(cycles)])
    all_cycles = np.zeros(size)
    MaxSize = min(cat_cyc.size,all_cycles.size)
    all_cycles[:MaxSize] = cat_cyc[:MaxSize]
    # XXX just extrapolate end..
    all_cycles[MaxSize:] = cat_cyc[-1]
    UnfoldObj = [RobTimeSepForceToIWT(o,ZFunc=(lambda: up)) for o in retracts]
    RefoldObj = [RobTimeSepForceToIWT(o,ZFunc=(lambda: down)) for o in reverse]
    NumericallyGetDeltaA(UnfoldObj,RefoldObj)
    # get the IWT
    Bins = [25,50,100,200,250,300,500,1000]
    for b in Bins:
        LandscapeObj =  InverseWeierstrass.\
                        FreeEnergyAtZeroForce(UnfoldObj,
                                              NumBins=b)
        fig = PlotUtilities.figure(figsize=(8,8))
        IWT_Util.ForceExtensionHistograms(Example.Separation,
                                          Example.Force,
                                          nBins=b)
        PlotUtilities.savefig(fig,OutBase + "0_{:d}hist.pdf".format(b))
        fig = PlotUtilities.figure(figsize=(12,12))
        IWT_Util.EnergyLandscapePlot(LandscapeObj,FOneHalf=16e-12,
                                     ZoomBoundsMeters=[56e-9,68e-9],
                                     stiffness_pN_per_nm=80,
                                     NumPointsAround=int(b/10))
        PlotUtilities.savefig(fig,OutBase + "1_{:d}IWT.pdf".format(b))
    
if __name__ == "__main__":
    run()
