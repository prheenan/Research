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
sys.path.append("../../../")
sys.path.append("../../../../../../../../../")
from LandscapeGeneration import LandscapeUtil
from Processing import Util
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,\
    FEC_Plot

from scipy.interpolate import LSQUnivariateSpline


def _get_region_demarcations(spline,time,first_N_peaks):
    """
    Gets the unfolding/refolding style demarcations
    
    Args:
        spline:
        time: the time bases to feed to the spline and its derivatives
        first_N_peaks: by default, we just want the first n peaks (the others
        are the actual unfolding peaks)
        
    Returns:
        tuple of <unfolding starts (size N+1), unfolding ends (size N)>. 
        There is one extra unfolding end to demarcate the end 
    """
    spline_t = spline(time)
    d_spline_dt = spline.derivative(n=1)(time)
    abs_d_spline_dt = np.abs(d_spline_dt)
    d2_spline_dt2 = spline.derivative(n=2)(time)
    negative_d2 = (d2_spline_dt2 < 0)
    zero_d1 = (abs_d_spline_dt < np.percentile(abs_d_spline_dt,3))
    where_turning_point = np.where(negative_d2 & zero_d1)[0]
    assert where_turning_point.size > 0 , "Couldn't find local maxima..."
    # get the centers of each turning point
    turning_point_idx = np.where(np.diff(where_turning_point) > 1)[0]
    # turning_point_idx is the list of indices where we jump from
    # peak a to b.
    start_idx = [0] + list(turning_point_idx+1)
    end_idx = list(turning_point_idx+1) + [-1]
    regions = [list(where_turning_point[i:f])
               for i,f in zip(start_idx,end_idx)]
    # make sure all the indices are there; have to add the last
    regions[-1] += [regions[-1][-1] + 1]
    assert (np.concatenate(regions) == where_turning_point).all() , \
        "Didn't use  them all"
    # determine where the derivative is locally minimized -- this is the
    # 'top' of the unfolding ramp, where we start to refold
    idx_refold_top = [r[0] + np.argmin(abs_d_spline_dt[r]) for r in regions]
    # take the first N peaks (only looking at the forward and reverse
    # experiments)
    assert len(idx_refold_top) >= first_N_peaks , "Didn't find enough peaks..."
    idx_refold_top = idx_refold_top[:first_N_peaks]
    # determine, in between those regions, where the spline is minimized
    # (this is the 'bottom' of the unfolding ramp
    idx_refold_bottom =  [i + np.argmin(spline_t[i:f])
                          for i,f in \
                          zip(idx_refold_top[:-1],idx_refold_top[1:])]
    # right now, the data looks like
    #  / \  / \  / \
    # but we are missing the first and the last points, '.' below: 
    # ./ \  / \  / \.
    # get the differences between the indices
    diff = lambda _i,_j : spline_t[_j] - spline_t[_i]
    diffs_fwd = [diff(i,j) for i, j in
                 zip(idx_refold_top[:-1], idx_refold_bottom)]
    diffs_rev = [diff(i,j) for i, j in
                 zip(idx_refold_top[1:], idx_refold_bottom)]
    all_diffs = np.concatenate([np.abs(diffs_fwd), np.abs(diffs_rev)])
    target_diff = min(all_diffs)
    # tack on the start and end of the experiment (both 'bottom' indices)
    # ie, the '.' above 
    offset_z = spline_t[idx_refold_top[0]]
    last_top = idx_refold_top[-1]
    start = np.where(spline_t[:last_top] >= offset_z - target_diff)[0][0]
    end = last_top + \
        np.where(spline_t[last_top:] <= offset_z - target_diff)[0][0]
    # add the two bottom indices
    idx_refold_bottom = [start] + idx_refold_bottom + [end]
    return idx_refold_bottom,idx_refold_top

def get_iwt_regions(input_dir,n_desired_resolution_m = 3e-9,N_peaks = 5):
    """
    Returns the *un-masked* regions (ie: not all the same length or DeltaZ)
    of input_dir. 
    """
    examples = CheckpointUtilities.lazy_multi_load(input_dir)
    for tmp in examples:
        e = tmp._slice(slice(0,None,1))
        z_raw = e.ZSnsr
        time = e.Time
        n_bins = int(np.ceil((max(z_raw)-min(z_raw))/n_desired_resolution_m))
        bins = np.linspace(min(time),max(time),endpoint=True,num=n_bins)
        spline = LSQUnivariateSpline(x=time,y=z_raw,t=bins[1:-1],k=3)
        spline_t = spline(time)

        # determine where all the events start and end 
        idx_refold_bottom,idx_refold_top = \
            _get_region_demarcations(spline,time,first_N_peaks=N_peaks)
        slice_experiment_region = slice(idx_refold_bottom[0],
                                        idx_refold_bottom[-1],None)
        # get the points *after* each turning point where we are
        # slice into the original data
        data = tmp._slice(slice_experiment_region)
        # determine the regions we care about
        to_ret = LandscapeUtil.RefoldingInfo(data,
                                             idx_start=idx_refold_bottom,
                                             idx_end=idx_refold_top,
                                             spline=spline)                    
        yield to_ret                                            

def run():
    input_dir = Util.cache_aligned("../../../")
    cache_dir = LandscapeUtil.cache_landscape_regions("../../../")
    load_f = lambda: get_iwt_regions(input_dir)
    e = CheckpointUtilities.multi_load(cache_dir=cache_dir,
                                       load_func=load_f,
                                       force=True,limit=None,
                                       name_func=FEC_Util.fec_name_func)
    objs = [tmp.data for tmp in e]
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(objs)
    PlotUtilities.savefig(fig,"./out")

if __name__ == "__main__":
    run()
