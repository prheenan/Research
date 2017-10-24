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


def get_aligned_regions(input_dir,v=50e-9):
    """
    Returns the 'aligned' iwt regions 
    """
    examples = CheckpointUtilities.lazy_multi_load(input_dir)
    spline_t = [e.spline(e.data.Time) for e in examples]
    # get the bounds for the region
    max_of_min = max([e._z_region_max_min for e in examples])
    # XXX why is this not min(max( ? 
    min_of_max = min([e._z_region_min_max for e in examples])
    len_expected = None
    key = examples[0]
    n_expected_pairs = len(key._idx_pairs)
    key_t = key.data.Time
    dt = key_t[1]-key_t[0]
    dz = abs(max_of_min-min_of_max)
    dn = int(np.ceil(dz/(v*dt)))
    for i,example_info in enumerate(examples):
        original = example_info.data._slice(slice(0,None,1))
        z = spline_t[i]
        t = example_info.data.Time
        slice_idx = example_info._idx_pairs
        z_schedule = np.zeros(t.size)
        unfolding = [original._slice(slice(r.start,r.start+dn,1))
                     for r in slice_idx[::2]]
        refolding = [original._slice(slice(r.start-dn,r.start,1))
                     for r in (slice_idx[2::2] + [slice(-1,0,1)])]
        for u,r in zip(unfolding,refolding):
            u.Velocity = v
            r.Velocity = -v
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
        # make sure the sizes match all previous
        if (len_expected is None):
            len_expected = key_length
        else:
            assert len_expected == key_length
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
