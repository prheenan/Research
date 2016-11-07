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

def RobTimeSepForceToIWT(o,ZFunc):
    # spring constant should be in N/m
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
    UnfoldObj = [RobTimeSepForceToIWT(o,ZFunc=(lambda: up))
                 for o in retracts[:1]]
    RefoldObj = [RobTimeSepForceToIWT(o,ZFunc=(lambda: down))
                 for o in reverse[:1]]
    # get the IWT
    Bins = [50,75,100,150,200,500,1000]
    kT = 4.1e-21
    for b in Bins:
        LandscapeObj =  InverseWeierstrass.\
            FreeEnergyAtZeroForce(UnfoldObj,NumBins=b,RefoldingObjs=RefoldObj)
        LandscapeObjFwd =  InverseWeierstrass.\
            FreeEnergyAtZeroForce(UnfoldObj,NumBins=b)
        # make a 2-D histogram of everything
        fig = PlotUtilities.figure(figsize=(8,8))
        ext_nm = Example.Separation*1e9
        IWT_Util.ForceExtensionHistograms(ext_nm,
                                          Example.Force*1e12,
                                          AddAverage=False,
                                          nBins=b)
        PlotUtilities.savefig(fig,OutBase + "0_{:d}hist.pdf".format(b))
        # get the distance to the transition state etc
        bounds_folded_nm = [56,61]
        bounds_unfolded_nm = [64,69]
        bounds_transition_nm = [61.75,63.5]
        FOneHalf_N = 15e-12
        # get the energy landscapes in kT
        Landscape_kT =  LandscapeObj.EnergyLandscape/kT
        Tilted_kT = Landscape_kT - \
                    (LandscapeObj.Extensions*FOneHalf_N)/kT
        landscape_ext_nm = LandscapeObj.Extensions * 1e9
        # fit a second order to the tilted one (easy to find transition states)
        pred_fold_x,pred_fold_y,coeffs_fold = \
            FitToRegion(landscape_ext_nm,Tilted_kT,bounds_folded_nm)
        pred_tx_x,pred_tx_y,coeffs_tx = \
            FitToRegion(landscape_ext_nm,Tilted_kT,bounds_transition_nm)
        pred_unfold_x,pred_unfold_y,coeffs_unfold = \
            FitToRegion(landscape_ext_nm,Tilted_kT,bounds_unfolded_nm)
        # get the offsets for the SHO
        x0_tx = ExtensionOffsetFromCoeffs(coeffs_tx)
        x0_fold = ExtensionOffsetFromCoeffs(coeffs_fold)
        x0_unfold = ExtensionOffsetFromCoeffs(coeffs_unfold)
        # DeltaX dagger is the distance between the transition and
        # folded state. See after equation 5:
        """
        Dudko, O. K., Hummer, G. & Szabo, A. 
        Theory, analysis, and interpretation of single-molecule force 
        spectroscopy experiments. PNAS 105, 1575515760 (2008).
        """
        DeltaXDagger = x0_tx-x0_fold
        print("DeltaX_Dagger is {:.1f}nm".format(DeltaXDagger))
        """
        DeltaG_dagger is (ibid)
        "the apparent free-energy of activation in the absence of an external
        force." (ie: the height of the energy barrier without force)
        """
        ext_idx = np.argmin(np.abs(landscape_ext_nm - x0_tx))
        DeltaGDagger = Landscape_kT[ext_idx]
        print("DeltaG_Dagger is {:.1f}kT".format(DeltaGDagger))
        fig = PlotUtilities.figure(figsize=(12,12))
        plt.subplot(2,1,1)
        plt.plot(landscape_ext_nm,Landscape_kT)
        plt.axvline(landscape_ext_nm[ext_idx],linewidth=4,color='g',
                    linestyle='--',
            label=(r"$\Delta x^{\dagger}$=" +
                   "{:.1f}nm".format(DeltaXDagger) ))
        plt.axhline(DeltaGDagger,linewidth=4,color='b',
                    linestyle='--',
                    label=(r"$\Delta G^{\dagger}}$=" +
                           "{:.1f}kT".format(DeltaGDagger) ))
        plt.ylim([0.5,max(Landscape_kT)*1.05])
        PlotUtilities.lazyLabel("","Landscape at F=0","",frameon=True)
        plt.subplot(2,1,2)
        Offset = np.percentile(Tilted_kT,10)
        OffsetTilted = Tilted_kT-Offset
        plt.plot(landscape_ext_nm,OffsetTilted,color='b',alpha=0.7)
        plt.axvline(landscape_ext_nm[ext_idx],linewidth=4,color='g',
                    linestyle='--',
            label=(r"$\Delta x^{\dagger}$=" +
                   "{:.1f}nm".format(DeltaXDagger) ))
        plt.plot(pred_fold_x,np.polyval(coeffs_fold,pred_fold_x)-Offset,
                 linestyle='--',color='r',linewidth=4,
                 label="Folded State at {:.1f}nm".format(x0_fold))
        plt.plot(pred_tx_x,np.polyval(coeffs_tx,pred_tx_x)-Offset,
                 linestyle='--',color='g',linewidth=4,
                 label="Transition State at {:.1f}nm".format(x0_tx))
        plt.plot(pred_unfold_x,np.polyval(coeffs_unfold,pred_unfold_x)-Offset,
                 linestyle='--',color='r',linewidth=4,
                 label="Unfolding State at {:.1f}nm".format(x0_unfold))
        plt.ylim(-0.5,max(OffsetTilted)*1.5)
        PlotUtilities.lazyLabel("Extension [nm]","Landscape at F1/2","",
                                frameon=True)
        PlotUtilities.savefig(fig,OutBase + "1_{:d}IWT.pdf".format(b))

def ExtensionOffsetFromCoeffs(coeffs):
    return -coeffs[1]/(2*coeffs[0])
        
def FitToRegion(x,y,bounds_x):
    GoodIdx = np.where( ( x > min(bounds_x)) &
                        ( x < max(bounds_x)) )
    pred_x = x[GoodIdx]
    pred_y = y[GoodIdx]
    coeffs = np.polyfit(x=pred_x,y=pred_y,deg=2)
    return pred_x,pred_y,coeffs

if __name__ == "__main__":
    run()
