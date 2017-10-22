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


def get_aligned_regions(input_dir):
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
    for i,example_info in enumerate(examples):
        original = example_info.data._slice(slice(0,None,1))
        z = spline_t[i]
        condition_between = ((z >= max_of_min) & \
                             (z <= min_of_max))
        slice_idx = example_info._idx_pairs
        idx_array = np.arange(0,z.size)
        regions = [np.where(condition_between & \
                            (idx_array >= r.start) & 
                            (idx_array <= r.stop))
                   for r in slice_idx]
        len_v = [len(r) for r in regions]
        # make sure there is at least one region
        assert (len_v > 0)
        # make sure all regions have the same length
        key_length = len_v[0]
        assert (len_v == [key_length for _ in len_v])
        # make sure these regions match all previous regions
        if (len_expected is None):
            len_expected = key_length
        else:
            assert len_expected == key_length
        to_ret = copy.deepcopy(example_info)
        # make a 
        data = [original._slice(r) for r in regions]
        unfolding = data[::2]
        refolding = data[1::2]
        assert len(unfolding) == len(refolding)
        assert (len(unfolding) + len(refolding)) == len(slice_idx)
        to_ret = LandscapeUtil.UnfoldingRefolding(unfolding,
                                                  refolding,
                                                  info=to_ret)
        yield to_ret

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
