# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.WormLikeChain.Python.Code import WLC_Utils
from scipy.interpolate import interp2d
from scipy.optimize import minimize_scalar

def optimize(L0,ext,f,cutoff=0.975,kbT=4.1e-21,Lp=0.33e-9):
    if (ext > L0 * cutoff):
        return f
    expected = WLC_Utils.WlcNonExtensible(ext=ext,kbT=kbT,Lp=Lp,L0=L0)
    return abs(f - expected)

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    _,data = FEC_Util.read_and_cache_pxp(directory="./",force=False)
    # offset and flip everything. All units are SI 
    x_offset = -4.62e-7
    f_offset = -1e-12
    max_ext = 60e-9
    # make a grid to convert to contour length space.
    ext = np.linspace(0,max_ext,num=max_ext * 1e9)
    min_L0 = 20e-9
    max_L0 = 70e-9
    fig = PlotUtilities.figure()
    for d in data:
        d.Force *= -1
        d.Force -= f_offset
        d.Separation -= x_offset
        plt.plot(d.Separation*1e9,d.Force*1e12,linewidth=0.5)
    # plot some wlcs...
    for i,L in enumerate(np.linspace(min_L0,max_L0,10)):
        f = WLC_Utils.WlcNonExtensible(ext=ext,kbT=4.1e-21,Lp=0.3e-9,L0=L)
        good_idx = np.where(ext < 0.9 * L)
        label = "wlc" if (i == 0) else ""
        plt.plot(ext[good_idx]*1e9,f[good_idx]*1e12,linewidth=0.75,alpha=0.5,
                 label=label)
    PlotUtilities.lazyLabel("Separation (nm)","Force (pN)","")
    PlotUtilities.savefig(fig,"fec.png")
    # look at the first data set. this is the first retraction,
    # which is *not* the actual flickering. the actual flickering goes to
    # zero a bunch, which messes up how we determine L0(F,q)
    key = data[0]
    idx_i = np.where(key.Separation > 0)[0][0]
    idx_f = np.where(key.Separation > max_L0)[0]
    if (idx_f.size > 0):
        idx_f = idx_f[0]
    else:
        idx_f = key.Separation.size
    # only look at the first N points; just to save time, since this is a demo
    idx_f = idx_i + min(idx_f-idx_i,int(1e4))
    m_slice = slice(idx_i,idx_f,1)
    # use a bounded root finder; for each pair of x,F, determine the closest
    # wlc curve, to an absolute tolerance of tol_m
    tol_m = 1e-12
    min_options = dict(bounds=(min_L0,max_L0),method="bounded",
                       options=dict(maxiter=int(1e3),disp=True,xatol=tol_m))
    L0_res = [minimize_scalar(fun=(lambda tmp: optimize(tmp,ext=ext,f=f)),
                              **min_options)
              for ext,f in zip(key.Separation[m_slice],key.Force[m_slice])]
    L0_opt = np.array([ l.x for l in L0_res])
    fig = PlotUtilities.figure((3.5,5))
    ax = plt.subplot(2,1,1)
    plt.plot(key.Time[m_slice],key.Force[m_slice]*1e12)
    PlotUtilities.lazyLabel("","Force (pN)","")
    PlotUtilities.no_x_label(ax)
    plt.subplot(2,1,2)
    plt.plot(key.Time[m_slice],L0_opt*1e9)
    PlotUtilities.lazyLabel("Time (s)","L$_0$ (nm)","")
    PlotUtilities.savefig(fig,"L0.png")


if __name__ == "__main__":
    run()
