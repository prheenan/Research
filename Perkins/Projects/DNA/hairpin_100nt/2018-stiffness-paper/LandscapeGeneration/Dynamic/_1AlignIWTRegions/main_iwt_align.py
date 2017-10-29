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

from scipy.interpolate import LSQUnivariateSpline


def get_unfolding_and_refolding(example_info,n_expected_pairs,dn,v):
    t = example_info.data.Time
    z_spline = example_info.spline
    z = z_spline(t)
    original = example_info.data._slice(slice(0,None,1))
    slice_idx = example_info._idx_pairs
    z_schedule = np.zeros(t.size)
    max_of_min = example_info._z_region_max_min
    z0 = max_of_min
    idx_starts = [r.start + np.where(z[r] >= z0)[0][0] for r in slice_idx]
    idx_start_unfold = idx_starts[::2]
    idx_end_refold = idx_starts[2::2] + [-1]
    unfolding = [original._slice(slice(r,r+dn,1))
                 for r in idx_start_unfold]
    refolding = [original._slice(slice(r-dn,r,1))
                 for r in idx_end_refold]
    for u,r in zip(unfolding,refolding):
        u_t_tmp = u.Time
        z = (u_t_tmp - u_t_tmp[0]) * v + z0
        u.Velocity = v
        u.ZSnsr = z.copy()
        r.Velocity = -v
        r.ZSnsr = z[::-1].copy()
        
    xlim = [min(original.ZSnsr),max(original.ZSnsr)]
    ylim = [min(original.Force),max(original.Force)]
    plt.subplot(2,1,1)
    plt.plot(original.ZSnsr,original.Force)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.subplot(2,1,2)
    for u,r in zip(unfolding,refolding):
        plt.plot(u.ZSnsr,u.Force)
        plt.plot(r.ZSnsr,r.Force)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()
    # make sure we have the same number of slices
    n_total_pairs = len(unfolding) + len(refolding)
    assert len(unfolding) == len(refolding)
    assert n_total_pairs == len(slice_idx)
    # make sure the number of pairs matches the expected number of pairs.
    assert n_total_pairs == n_expected_pairs
    # make sure we have the same sizes for all sizes
    _check_sizes_consistent(unfolding)
    _check_sizes_consistent(refolding)
    # check that the two lists are consistent
    key_length = unfolding[0].Force.size
    assert key_length == refolding[0].Force.size
    return unfolding,refolding

def get_n_points(key,v):
    n_expected_pairs = len(key._idx_pairs)
    key_t = key.data.Time
    dt = key_t[1]-key_t[0]
    max_of_min = key._z_region_max_min
    min_of_max = key._z_region_min_max
    dz = abs(max_of_min-min_of_max)
    dn = int(np.ceil(dz/(v*dt)))
    return dn


def get_aligned_regions(input_dir,v=50e-9):
    """
    Returns the 'aligned' iwt regions 
    """
    examples = CheckpointUtilities.lazy_multi_load(input_dir)
    n_expected_pairs = 10
    for i,example_info in enumerate(examples):
        dn = get_n_points(example_info,v)
        unfolding,refolding = \
            get_unfolding_and_refolding(example_info,n_expected_pairs,dn,v)
        # make sure the sizes match all previous
        to_ret = copy.deepcopy(example_info)
        to_ret = LandscapeUtil.UnfoldingRefolding(unfolding,
                                                  refolding,
                                                  info=to_ret)
        yield to_ret

def _check_sizes_consistent(list_v):
    key = list_v[0].Force.size
    sizes = [u.Force.size for u in list_v]
    err_string = "Actual sizes ({:s}) don't match the key ({:d})".\
                 format(sizes,key)
    assert sizes == [key for _ in list_v] , err_string


def run():
    input_dir =  LandscapeUtil.cache_landscape_regions("../../../")
    cache_dir = LandscapeUtil.cache_landscape_aligned_regions("../../../")
    load_f = lambda: get_aligned_regions(input_dir)
    name_func = lambda i,e: FEC_Util.fec_name_func(i,e) + "{:d}".format(i)
    e = CheckpointUtilities.multi_load(cache_dir=cache_dir,
                                       load_func=load_f,
                                       force=True,limit=None,
                                       name_func=FEC_Util.fec_name_func)
    fecs_unfold = [d for tmp in e for d in tmp.unfolding]
    fecs_refold = [d for tmp in e for d in tmp.refolding]
    for u,r in zip(fecs_unfold,fecs_refold):
        unfold_size = u.Force.size
        refold_size = r.Force.size
        error_str = "Sizes ({:d}/{:d}) don't match".format(unfold_size,
                                                           refold_size)
        assert unfold_size == refold_size , error_str
    all_sep = np.concatenate([d.Separation 
                              for d in (fecs_unfold + fecs_refold)])
    min_x = np.min(all_sep)
    max_x = np.max(all_sep)
    xlim = [min_x * 1e9, max_x * 1e9]
    fig = PlotUtilities.figure((4,8))
    # plot the heatmaps...
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
    fig = PlotUtilities.figure()
    for u,r in zip(fecs_unfold,fecs_refold):
        plt.plot(u.ZSnsr*1e9,u.Force*1e12)
        plt.plot(r.ZSnsr*1e9,r.Force*1e12)
    PlotUtilities.lazyLabel("Z (nm)","Force (pN)","")
    PlotUtilities.savefig(fig,"./F_vs_z")


if __name__ == "__main__":
    run()
