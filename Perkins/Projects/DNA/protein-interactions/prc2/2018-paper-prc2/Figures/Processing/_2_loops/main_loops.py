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
from Util import IoUtil
from Research.Perkins.AnalysisUtil.Images import PolymerTracing
from GeneralUtil.python import CheckpointUtilities
from scipy.interpolate import splev,BSpline
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

def detect_loops_in_trace(wlc):
    x,y = wlc._x_raw,wlc._y_raw
    u,tck,spline,deriv = PolymerTracing._u_tck_spline_and_derivative(x,y)
    # get the 
    n_points_interp = u.size/len(x)
    assert abs(n_points_interp - int(n_points_interp)) < 1e-6
    # POST: n_points_interp is an integer
    tolerance = 0.1
    n_points_interp = int(np.ceil(3/tolerance))
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
    pairs_of_splines = set()
    for i,(min_idx,s1) in enumerate(zip(min_dist_idx,spline_slices)):
        # make a minimize function for this spline
        # only ask for a solution within <xatol> pixels
        n_points = max(2,int(2*np.ceil(1 / tolerance)))
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
            pairs_of_splines.add(pair)
    plt.close()
    plt.plot(x,y,'b.-')
    for s in spline_slices:
        plt.plot(*s)
    debug_idx = [4,14,23,35]
    for i_tmp in debug_idx:
        plt.plot(*spline_slices[i_tmp],color='r',linewidth=3)
        plt.plot(*(spline_slices[i_tmp][:,0]),color='r',marker='o')
        plt.plot(*spline_slices[min_dist_idx[i_tmp]],color='r',linewidth=3)
    for p in parameter_values_where_crossover:
        plt.plot(*PolymerTracing.evaluate_spline(p, tck=tck),
                 color='k',marker='x')
    plt.show()
    print("...")

def run(in_dir):
    """
    Args:
        in_dir: the input directory to operate on.  
    """
    input_dir =  IoUtil.data_dir(in_dir)
    cache_dir = IoUtil._traces_dir(in_dir)
    output_dir = IoUtil._ensemble_dir(in_dir)
    # just read in from the cache...
    objs_all = IoUtil.read_images(input_dir,cache_dir=cache_dir,force=False,
                                  limit=1)
    tmp = objs_all[0]
    wlc = tmp.worm_objects[1]
    detect_loops_in_trace(wlc)
    idx_debug = 15
    x_raw,y_raw = wlc._x_raw,wlc._y_raw
    x,y = wlc.inf.x_y_abs
    
    
if __name__ == "__main__":
    run(IoUtil.get_directory_command_line())
