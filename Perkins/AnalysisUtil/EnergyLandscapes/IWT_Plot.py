# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from GeneralUtil.python import PlotUtilities
import IWT_Util
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass


def TomPlot(LandscapeObj,OutBase,UnfoldObj,RefoldObj,idx,bounds):
    # get a forward and reverse
    ToX = lambda x: x * 1e9
    ToForceY = lambda y: y * 1e12
    fig = PlotUtilities.figure(figsize=(8,4))
    plt.subplot(1,2,1)
    SubplotArgs = dict(alpha=0.4,linewidth=0.5)
    FilterN = 500
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
    Obj =  IWT_Util.TiltedLandscape(LandscapeObj,bounds)
    plt.plot(Obj.landscape_ext_nm,Obj.OffsetTilted_kT)
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
    np.savetxt(X=np.c_[Obj.landscape_ext_nm,Obj.OffsetTilted_kT],
               fname=OutBase+"Landscape"+ext,**common)
    
def plot_single_landscape(LandscapeObj,bounds=None,min_landscape_kT=None,
                          max_landscape_kT=None):
    """
    Plots a detailed energy landscape, and saves

    Args:
        LandscapeObj: energy landscape object (untilted)
        bounds: Where to 'zoom' in the plot, Iwt_Util.BoundsObj instance
        Bins: how many bins to use in the energy landscape plots
        <min/max>_landscape_kT: bounds on the landscape
    Returns:
        nothing
    """                          
    if (bounds is None):
        all = [0,np.inf]
        bounds = IWT_Util.BoundsObj(all,all,all,all)
    Obj =  IWT_Util.TiltedLandscape(LandscapeObj,
                                    bounds)
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
    plt.plot(Obj.landscape_ext_nm,Obj.OffsetTilted_kT,color='b',alpha=0.7)
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
    if (max_landscape_kT is None):
        max_landscape_kT = max(Obj.OffsetTilted_kT)*1.5
    if (min_landscape_kT is None):
        min_landscape_kT = np.percentile(Obj.OffsetTilted_kT,5)-2
    plt.ylim( min_landscape_kT,max_landscape_kT)
    PlotUtilities.lazyLabel("Extension [nm]","Landscape at F1/2","",
                            frameon=True)
                            
def InTheWeedsPlot(OutBase,UnfoldObj,bounds=None,RefoldObj=[],Example=None,
                   Bins=[50,75,100,150,200,500,1000],**kwargs):
    """
    Plots a detailed energy landscape, and saves

    Args:
        OutBase: where to start the save
        UnfoldObj: unfolding objects
        bounds: Where to 'zoom' in the plot, Iwt_Util.BoundsObj instance
        RefoldObj: refolding objects
        Bins: how many bins to use in the energy landscape plots
        <min/max>_landscape_kT: bounds on the landscape
    Returns:
        nothing
    """
    # get the IWT
    kT = 4.1e-21
    for b in Bins:
        LandscapeObj =  InverseWeierstrass.\
            FreeEnergyAtZeroForce(UnfoldObj,NumBins=b,RefoldingObjs=RefoldObj)
        # make a 2-D histogram of everything
        if (Example is not None):
            fig = PlotUtilities.figure(figsize=(8,8))
            ext_nm = Example.Separation*1e9
            IWT_Util.ForceExtensionHistograms(ext_nm,
                                              Example.Force*1e12,
                                              AddAverage=False,
                                              nBins=b)
            PlotUtilities.savefig(fig,OutBase + "0_{:d}hist.pdf".format(b))
        # get the distance to the transition state etc

        print("DeltaG_Dagger is {:.1f}kT".format(Obj.DeltaGDagger))
        fig = PlotUtilities.figure(figsize=(12,12))
        plot_single_landscape(LandscapeObj,bounds=bounds,**kwargs)
        PlotUtilities.savefig(fig,OutBase + "1_{:d}IWT.pdf".format(b))

