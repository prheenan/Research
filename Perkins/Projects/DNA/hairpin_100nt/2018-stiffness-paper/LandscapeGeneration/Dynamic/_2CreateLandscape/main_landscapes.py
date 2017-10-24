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


def get_aligned_regions(input_dir):
    """
    Returns the 'aligned' iwt regions 
    """

def run():
    input_dir =  LandscapeUtil.cache_landscape_aligned_regions("../../../")
    cache_dir = LandscapeUtil.cache_landscapes("../../../")
    e = CheckpointUtilities.lazy_multi_load(input_dir)
    fecs_unfold = [WeierstrassUtil.ToIWTObject(d) for tmp in e 
                   for d in tmp.unfolding]
    fecs_refold = [WeierstrassUtil.ToIWTObject(d) for tmp in e 
                   for d in tmp.refolding]
    z0 = 20e-9
    key = fecs_unfold[0]
    t = key.Time
    zf = z0 + (t[-1]-t[0]) * key.Velocity
    for un,re in zip(fecs_unfold,fecs_refold):
        un.SetOffsetAndVelocity(Offset=z0,Velocity=un.Velocity)
        re.SetOffsetAndVelocity(Offset=zf,Velocity=re.Velocity)
    l = InverseWeierstrass.\
        free_energy_inverse_weierstrass(fecs_unfold,fecs_refold)
    g_rel_z_offset_m = min(l.q)
    stiffness_rel = 2e-12/3e-9
    g_relative = (stiffness_rel/2) * (l.q - g_rel_z_offset_m)**2 
    G0 = l.G_0 - min(l.G_0)
    plt.subplot(2,1,1)
    plt.plot(l.q,G0)
    plt.plot(l.q,g_relative)
    plt.subplot(2,1,2)
    plt.plot(l.q,G0-g_relative)
    plt.show()
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
