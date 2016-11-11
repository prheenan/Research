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
# file-wide bounds on our state extensions in nm
bounds_folded_nm = [56,61]
bounds_unfolded_nm = [64,69]
bounds_transition_nm = [61.75,63.5]
FOneHalf_N = 14.2e-12



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

class TiltedLandscape:
    def __init__(self,landscape,bounds_folded_nm,bounds_transition_nm,
                 bounds_unfolded_nm,f_one_half,kT=4.1e-21):
        self.Landscape_kT =  landscape.EnergyLandscape/kT
        self.Tilted_kT = self.Landscape_kT - \
                         (landscape.Extensions*f_one_half)/kT
        landscape_ext_nm = landscape.Extensions * 1e9
        self.pred_fold_x,self.pred_fold_y,self.coeffs_fold = \
            FitToRegion(landscape_ext_nm,self.Tilted_kT,bounds_folded_nm)
        self.pred_tx_x,self.pred_tx_y,self.coeffs_tx = \
                FitToRegion(landscape_ext_nm,self.Tilted_kT,
                            bounds_transition_nm)
        self.pred_unfold_x,self.pred_unfold_y,self.coeffs_unfold = \
                FitToRegion(landscape_ext_nm,self.Tilted_kT,bounds_unfolded_nm)
        # get the energy landscapes in kT
        # fit a second order to the tilted one (easy to find transition states)
        # get the offsets for the SHO
        self.x0_tx = ExtensionOffsetFromCoeffs(self.coeffs_tx)
        self.x0_fold = ExtensionOffsetFromCoeffs(self.coeffs_fold)
        self.x0_unfold = ExtensionOffsetFromCoeffs(self.coeffs_unfold)
        # DeltaX dagger ais the distance between the transition and
        # folded state. See after equation 5:

        """
        Dudko, O. K., Hummer, G. & Szabo, A. 
        Theory, analysis, and interpretation of single-molecule force 
        spectroscopy experiments. PNAS 105, 1575515760 (2008).
        """
        self.DeltaXDagger = self.x0_tx-self.x0_fold
        """
        DeltaG_dagger is (ibid)
        "the apparent free-energy of activation in the absence of an external
        force." (ie: the height of the energy barrier without force)
        """
        self.ext_idx = np.argmin(np.abs(landscape_ext_nm - self.x0_tx))
        self.MinG = min(self.Landscape_kT)
        self.DeltaGDagger = self.Landscape_kT[self.ext_idx]-self.MinG
        self.landscape_ext_nm = landscape_ext_nm
        self.Offset = np.percentile(self.Tilted_kT,10)
        self.OffsetTilted = self.Tilted_kT-self.Offset


def TomPlot(OutBase,UnfoldObj,RefoldObj,Bins,idx):
    LandscapeObj =  InverseWeierstrass.\
        FreeEnergyAtZeroForce(UnfoldObj,NumBins=Bins,RefoldingObjs=RefoldObj)
    # get a forward and reverse
    ToX = lambda x: x * 1e9
    ToForceY = lambda y: y * 1e12
    fig = PlotUtilities.figure(figsize=(8,4))
    plt.subplot(1,2,1)
    SubplotArgs = dict(alpha=0.4,linewidth=0.5)
    FilterN = 250
    Unfold = FEC_Util.GetFilteredForce(UnfoldObj[idx],FilterN)
    Refold = FEC_Util.GetFilteredForce(RefoldObj[idx],FilterN)
    UnfoldX = ToX(Unfold.Extension)
    UnfoldY = ToForceY(Unfold.Force)
    FoldX = ToX(Refold.Extension)
    FoldY = ToForceY(Refold.Force)
    plt.plot(UnfoldX,UnfoldY,color='r',label="Unfolding",
             **SubplotArgs)
    plt.plot(FoldX,FoldY,color='b',label="Refolding",
             **SubplotArgs)
    fontdict = dict(fontsize=13)
    x_text_dict =  dict(x=60, y=22.5, s="2 nm", fontdict=fontdict,
                        withdash=False,
                        rotation="horizontal")
    y_text_dict =  dict(x=59, y=27, s="5 pN", fontdict=fontdict, withdash=False,
                        rotation="vertical")
    PlotUtilities.ScaleBar(x_kwargs=dict(x=[60,62],y=[24,24]),
                           y_kwargs=dict(x=[60,60],y=[25,30]),
                           text_x=x_text_dict,text_y=y_text_dict)
    PlotUtilities.legend(loc=[0.4,0.8],**fontdict)
    plt.subplot(1,2,2)
    Obj = TiltedLandscape(LandscapeObj,
                          bounds_folded_nm,
                          bounds_transition_nm,
                          bounds_unfolded_nm,
                          FOneHalf_N)
    plt.plot(Obj.landscape_ext_nm,Obj.OffsetTilted)
    plt.xlim([56,69])
    plt.ylim([-1,4])
    yoffset = 1
    x_text_dict =  dict(x=58.5, y=yoffset+1.5, s="2 nm", fontdict=fontdict,
                        withdash=False,rotation="horizontal")
    y_text_dict =  dict(x=57, y=yoffset+2.5, s=r"1 k$_\mathrm{b}$T",
                        fontdict=fontdict, withdash=False,
                        rotation="vertical")
    PlotUtilities.ScaleBar(x_kwargs=dict(x=[58,60],
                                         y=[yoffset+1.75,yoffset+1.75]),
                           y_kwargs=dict(x=[58,58],y=[yoffset+2,yoffset+3]),
                           text_x=x_text_dict,text_y=y_text_dict,
                           kill_axis=True)
    PlotUtilities.savefig(fig,OutBase + "TomMockup" + str(idx) + ".png",
                          subplots_adjust=dict(bottom=-0.1))
    # save out the data exactly as we want to plot it
    common = dict(delimiter=",")
    ext = str(idx) + ".txt"
    np.savetxt(X=np.c_[UnfoldX,UnfoldY],fname=OutBase+"Unfold" + ext,**common)
    np.savetxt(X=np.c_[FoldX,FoldY],fname=OutBase+"Fold"+ext,**common)
    np.savetxt(X=np.c_[Obj.landscape_ext_nm,Obj.OffsetTilted],
               fname=OutBase+"Landscape"+ext,**common)
    


def InTheWeedsPlot(OutBase,UnfoldObj,RefoldObj,Example,
                   Bins=[50,75,100,150,200,500,1000]):
    # get the IWT
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
        Obj = TiltedLandscape(LandscapeObj,
                              bounds_folded_nm,
                              bounds_transition_nm,
                              bounds_unfolded_nm,
                              FOneHalf_N)
        print("DeltaG_Dagger is {:.1f}kT".format(Obj.DeltaGDagger))
        fig = PlotUtilities.figure(figsize=(12,12))
        plt.subplot(2,1,1)
        plt.plot(Obj.landscape_ext_nm,Obj.Landscape_kT)
        plt.axvline(Obj.landscape_ext_nm[Obj.ext_idx],linewidth=4,color='g',
                    linestyle='--',
            label=(r"$\Delta x^{\ddag}$=" +
                   "{:.1f}nm".format(Obj.DeltaXDagger) ))
        plt.axhline(Obj.DeltaGDagger,linewidth=4,color='b',
                    linestyle='--',
                    label=(r"$\Delta G^{\ddag}}$=" +
                           "{:.1f}kT".format(Obj.DeltaGDagger) ))
        plt.ylim([-0.5,max(Obj.Landscape_kT)*1.05])
        PlotUtilities.lazyLabel("","Landscape at F=0","",frameon=True)
        plt.subplot(2,1,2)
        plt.plot(Obj.landscape_ext_nm,Obj.OffsetTilted,color='b',alpha=0.7)
        plt.axvline(Obj.landscape_ext_nm[Obj.ext_idx],linewidth=4,color='g',
                    linestyle='--',
            label=(r"$\Delta x^{\ddag}$=" +
                   "{:.1f}nm".format(Obj.DeltaXDagger) ))
        plt.plot(Obj.pred_fold_x,
                 np.polyval(Obj.coeffs_fold,Obj.pred_fold_x)-Obj.Offset,
                 linestyle='--',color='r',linewidth=4,
                 label="Folded State at {:.1f}nm".format(Obj.x0_fold))
        plt.plot(Obj.pred_tx_x,
                 np.polyval(Obj.coeffs_tx,Obj.pred_tx_x)-Obj.Offset,
                 linestyle='--',color='g',linewidth=4,
                 label="Transition State at {:.1f}nm".format(Obj.x0_tx))
        plt.plot(Obj.pred_unfold_x,
                 np.polyval(Obj.coeffs_unfold,Obj.pred_unfold_x)-Obj.Offset,
                 linestyle='--',color='r',linewidth=4,
                 label="Unfolding State at {:.1f}nm".format(Obj.x0_unfold))
        plt.ylim(-0.5,max(Obj.OffsetTilted)*1.5)
        PlotUtilities.lazyLabel("Extension [nm]","Landscape at F1/2","",
                                frameon=True)
        PlotUtilities.savefig(fig,OutBase + "1_{:d}IWT.pdf".format(b))
        
def ExtensionOffsetFromCoeffs(coeffs):
    """
    Gets the location of the center of the 'well' from the polynomial coeffs

    Args:
        coeffs: the polynomial coefficients, higherst first, from np.polyfit
    Returns:
        x0, from k/2 * (x-x0)**2, fit 
    """
    return -coeffs[1]/(2*coeffs[0])
        
def FitToRegion(x,y,bounds_x):
    """
    Fits to a bounded region

    Args:
        x,y: values to fit
        bounds_x: only fit within these bounds
    Returns:
        tuple of <predicted x, predicted y, coefficients>
    """
    GoodIdx = np.where( ( x >= min(bounds_x)) &
                        ( x <= max(bounds_x)) )
    pred_x = x[GoodIdx]
    pred_y = y[GoodIdx]
    coeffs = np.polyfit(x=pred_x,y=pred_y,deg=2)
    return pred_x,pred_y,coeffs

    
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
                 for o in retracts]
    RefoldObj = [RobTimeSepForceToIWT(o,ZFunc=(lambda: down))
                 for o in reverse]
    for idx in range(min(len(UnfoldObj),RefoldObj)):
        TomPlot(OutBase+ str(idx) + "_",UnfoldObj,RefoldObj,Bins=40,idx=idx)
    InTheWeedsPlot(OutBase ,UnfoldObj,RefoldObj,Example,
                   Bins = [50,75,100,150])


if __name__ == "__main__":
    run()
