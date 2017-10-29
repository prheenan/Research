# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,copy
sys.path.append("../../../")
sys.path.append("../../../../../../../../../")
from LandscapeGeneration import LandscapeUtil
from Processing import Util
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass,WeierstrassUtil

from scipy.interpolate import LSQUnivariateSpline


def convert_to_iwt(tmp):
    unfold = [WeierstrassUtil.ToIWTObject(d,Offset=d.ZSnsr[0])
              for d in tmp.unfolding]
    refold = [WeierstrassUtil.ToIWTObject(d,Offset=d.ZSnsr[0])
              for d in tmp.refolding]
    key = tmp.unfolding[0]
    z0 = key.ZSnsr[0]
    v = key.Velocity
    for u,r in zip(unfold,refold):
        t = u.Time
        zf = z0 + v * (t[-1]-t[0])
        u.SetOffsetAndVelocity(z0,+v)
        r.SetOffsetAndVelocity(zf,-v)
        u.Time -= min(u.Time)
        r.Time -= min(r.Time)
    return unfold,refold 

def get_landscapes(unfold,refold):
    l = InverseWeierstrass.free_energy_inverse_weierstrass(unfold)
    return l


def run():
    input_dir =  LandscapeUtil.cache_landscape_aligned_regions("../../../")
    cache_dir = LandscapeUtil.cache_landscapes("../../../")
    e = CheckpointUtilities.lazy_multi_load(input_dir)
    data = [convert_to_iwt(tmp) for tmp in e]
    landscapes = [get_landscapes(*d) for d in data]
    kw = dict(linewidth=0.5)
    for tmp in e:
        for u,r in zip(tmp.unfolding,tmp.refolding):
            plt.plot(u.ZSnsr,u.Force)
            plt.plot(r.ZSnsr,r.Force)
        plt.show()
    for i,l in enumerate(landscapes):
        plt.plot(l.q,l.G_0)
    plt.show()
    fecs_unfold = [WeierstrassUtil.ToIWTObject(d) for tmp in e 
                   for d in tmp.unfolding]
    fecs_refold = [WeierstrassUtil.ToIWTObject(d) for tmp in e 
                   for d in tmp.refolding]
    all_sep = np.concatenate([d.Separation 
                              for d in (fecs_unfold + fecs_refold)])
    min_x = np.min(all_sep)
    max_x = np.max(all_sep)
    xlim = [min_x * 1e9, max_x * 1e9]
    fig = PlotUtilities.figure((4,8))
    ax = plt.subplot(2,1,1)
    FEC_Plot.heat_map_fec(fecs_unfold)
    PlotUtilities.no_x_label(ax)
    PlotUtilities.xlabel("")
    plt.xlim(xlim)
    plt.subplot(2,1,2)
    FEC_Plot.heat_map_fec(fecs_refold)
    PlotUtilities.title("")
    plt.xlim(xlim)
    PlotUtilities.savefig(fig,"./out")


if __name__ == "__main__":
    run()
