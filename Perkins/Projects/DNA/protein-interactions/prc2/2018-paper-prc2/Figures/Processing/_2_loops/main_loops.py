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

sys.path.append("../../../../../../../../../../")
sys.path.append("../../")
from Util import IoUtil,AnalysisClasses
from Research.Perkins.AnalysisUtil.Images import PolymerTracing
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
from scipy.optimize import brute

def _root_func(i,s1,tck,max_n):
    try:
        u = i[0]
    except IndexError:
        u = i
    if ((i > max_n) or (u <= 1e-200)):
        return np.inf
    x_y_point = np.array(PolymerTracing.evaluate_spline(u, tck=tck))
    diff = (s1.T - x_y_point).T
    diff_abs_summed = np.sum(np.abs(diff),axis=0)
    to_ret = np.min(diff_abs_summed)
    return to_ret



def get_distances(x):
    """
    :param x: a list of spline objects representing the curve
    :return: distance[i][j] is the smallest distance from x[i] to x[j]
    """
    distances = np.zeros((len(x),len(x)))
    for i,s1 in enumerate(x):
        tmp = []
        for j,s2 in enumerate(x):
            if (j < i+2 and j > i-2):
                add = np.inf
            else:
                x1,x2 = s1[0,:],s2[0,:]
                y1,y2 = s1[1,:],s2[1,:]
                diff_matrix_x = PolymerTracing._difference_matrix(x1,x2)**2
                diff_matrix_y = PolymerTracing._difference_matrix(y1,y2)**2
                dist_sq = (diff_matrix_x + diff_matrix_y)
                add = np.min(np.sqrt(dist_sq))
            tmp.append(add)
        distances[i,:] = tmp
    return distances

def detect_loops_in_trace(wlc,upscale_interp=3,tolerance=0.1):
    x,y = wlc._x_raw,wlc._y_raw
    u,tck,spline,deriv = PolymerTracing._u_tck_spline_and_derivative(x,y)
    # get the 
    n_points_interp = u.size/len(x)
    assert abs(n_points_interp - int(n_points_interp)) < 1e-6
    # POST: n_points_interp is an integer
    n_points_interp = int(np.ceil(upscale_interp/tolerance))
    # u runs from 0 to n-1; get just the slice from i to i+1, for i in 0 to n-2
    # (since for N data points, there are only N-1 segments. Think N=2 -- just
    # one (N-1) line connecting two (N) points)
    get_u = lambda i : np.linspace(i,i+1,n_points_interp)
    spline_slices = [np.array(PolymerTracing.evaluate_spline(get_u(i),tck=tck))
                     for i,_ in enumerate(x[:-1])]
    n_spline_slices = len(spline_slices)
    distances = get_distances(spline_slices)
    # make sure each element compares to all (N-1) other elements
    np.testing.assert_allclose(distances.shape,\
                               (n_spline_slices,n_spline_slices))
    # min_dist_idx[i] is the closest idx to of spline_slices[i]
    min_dist_idx = [ np.nanargmin(d) for d in distances]
    # loop through each pair of splines, and determine the roots (closer than
    # <tolerance> of a pixel is considered an intersection).
    parameter_values_where_crossover = []
    pairs_of_splines = []
    for i,(min_idx,s1) in enumerate(zip(min_dist_idx,spline_slices)):
        # make a minimize function for this spline
        # only ask for a solution within <xatol> pixels
        n_points = max(2,int(upscale_interp*np.ceil(1 / tolerance)))
        step = 1/n_points
        bounds = [slice(min_idx,min_idx+1+step,step)]
        kw = dict(ranges=bounds,
                  Ns=n_points,
                  full_output=True)
        f_min_spline = lambda i_lambda: _root_func(i_lambda,s1,tck,
                                                   max_n=n_spline_slices)
        x0,fval,_,_ = brute(f_min_spline,**kw)
        # get the parameter difference between each..
        u_x = x0[0]
        # if the difference minimizing x and y is small enough, and
        # the minimization suceeded, and the minimization was small enough,
        # then we count this as a crossover
        value_func = fval
        is_small_value = (value_func < tolerance).all()
        is_in_range = (u_x <= min_idx + 1) and (u_x >= min_idx)
        pair = "".join("{:d}".format(d) for d in sorted([i, min_idx]))
        if (is_small_value and is_in_range):
            parameter_values_where_crossover.append(u_x)
            pairs_of_splines.append(pair)
    # get the start and end of each loop
    set_of_pairs = sorted(list(set(pairs_of_splines)))
    loop_starts,loop_ends = [],[]
    for s in set_of_pairs:
        # each look should stop and end.
        # make sure there are
        u_idx = [i for i,s_tmp in enumerate(pairs_of_splines) if s_tmp == s]
        u_loop = [parameter_values_where_crossover[i] for i in u_idx]
        # by definition, u is along the contour, so loop starts at min, ends
        # at max
        u0,u1 = min(u_loop),max(u_loop)
        loop_starts.append(u0)
        loop_ends.append(u1)
    return AnalysisClasses.Loops(u,tck,loop_starts,loop_ends)

def debug_plot(tmp_loop):
    spline_x_y = PolymerTracing.evaluate_spline(tmp_loop.u, tmp_loop.tck)
    plt.plot(*spline_x_y,label="Full DNA")
    for i,loop_bounds in enumerate(tmp_loop.get_loop_bounds()):
        u_min, u_max = loop_bounds[0], loop_bounds[1]
        num = (u_max - u_min) * 10
        if (u_max == u_min):
            print("uhh... umax=umin...")
            continue
        u = np.linspace(u_min, u_max, endpoint=True, num=num)
        spline_loop = PolymerTracing.evaluate_spline(u, tmp_loop.tck)
        label = "loop" if i == 0 else ""
        plt.plot(*spline_loop, linewidth=3,alpha=0.5,label=label)
    PlotUtilities.lazyLabel("x (pixels)","y (pixels)","")


def get_subset_loops(subset):
    loops = [detect_loops_in_trace(wlc) for wlc in subset]
    return loops

def run(in_dir):
    """
    Args:
        in_dir: the input directory to operate on.  
    """
    debug = False
    input_dir =  IoUtil.data_dir(in_dir)
    cache_dir = IoUtil._traces_dir(in_dir)
    output_dir = IoUtil._loop_dir(in_dir)
    # just read in from the cache...
    objs_all = IoUtil.read_images(input_dir,cache_dir=cache_dir,force=False)
    debug = True
    ensemble = AnalysisClasses.convert_to_ensemble(objs_all)
    # remove the multimers and unknown; we score them separately.
    ensemble.dna_only = get_subset_loops(ensemble.dna_only)
    ensemble.dna_plus_protein = get_subset_loops(ensemble.dna_plus_protein)
    ensemble.multimer = None
    ensemble.unknown = None
    # save out the ensemble
    CheckpointUtilities.lazy_save(output_dir + "loops.pkl",ensemble)
    if (debug):
        all_plottable = ensemble.dna_only + ensemble.dna_plus_protein
        for i,tmp_loop in enumerate(all_plottable):
            fig = PlotUtilities.figure()
            debug_plot(tmp_loop)
            PlotUtilities.savefig(fig,"./{:d}.png".format(i))

if __name__ == "__main__":
    run(IoUtil.get_directory_command_line())
