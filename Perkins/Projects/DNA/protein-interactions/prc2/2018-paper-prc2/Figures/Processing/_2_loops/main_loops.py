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
from scipy.optimize import root

def _root_func(i,s1,tck,max_n):
    if (i > max_n):
        return [np.inf,np.inf]
    diff = (s1.T - np.array(PolymerTracing.evaluate_spline(i[0], tck=tck))).T
    diff_summed = (np.abs(np.sum(diff,axis=0)))
    min_summed = np.argmin(diff_summed)
    return diff[:,min_summed]

def detect_loops_in_trace(wlc):
    x,y = wlc._x_raw,wlc._y_raw
    u,tck,spline,deriv = PolymerTracing._u_tck_spline_and_derivative(x,y)
    # get the 
    n_points_interp = u.size/len(x)
    assert abs(n_points_interp - int(n_points_interp)) < 1e-6
    # POST: n_points_interp is an integer
    tolerance = 0.25
    n_points_interp = int(10/tolerance)
    # u runs from 0 to n-1; get just the slice from i to i+1, for i in 0 to n-2
    # (since for N data points, there are only N-1 segments. Think N=2 -- just
    # one (N-1) line connecting two (N) points)
    get_u = lambda i : np.linspace(i,i+1,n_points_interp)
    spline_slices = [np.array(PolymerTracing.evaluate_spline(get_u(i),tck=tck))
                     for i,_ in enumerate(x[:-1])]
    n_spline_slices = len(spline_slices)
    # find the other segment which is closest in L1; anything crossing should 
    # be closest. Note that we sum over x and y. We ignore i == j, to prevent
    # the trivial case.
    distances = [ [np.sum(np.abs(s2-s1),axis=0) if (j >= i+2 or j <= i-2) else
                   [np.inf for _ in s1[0,:]]
                   for j,s1 in enumerate(spline_slices)] 
                  for i,s2 in enumerate(spline_slices)]
    # make sure each element compares to all (N-1) other elements
    np.testing.assert_allclose([len(d) for d in distances],n_spline_slices)
    # make sure each element is a list of the distances for each point
    np.testing.assert_allclose([len(s) for d in distances for s in d],
                               n_points_interp)
    # min_dist_idx[i] is the closest idx to of spline_slices[i]
    min_dist_idx = [ np.argmin([min(s) for s in d]) for d in distances]
    # loop through each pair of splines, and determine the roots (closer than
    # <tolerance> of a pixel is considered an intersection).
    parameter_values_where_crossover = []
    for i,(min_idx,s1) in enumerate(zip(min_dist_idx,spline_slices)):
        # make a minimize function for this spline
        # only ask for a solution within <xatol> pixels
        bounds = [min_idx,min_idx+1]
        kw = dict(x0 = min_idx + 1/2,
                  method="hybr",
                  options=dict(factor=1))
        f_min_spline = lambda i_lambda: _root_func(i_lambda,s1,tck,
                                                   max_n=n_spline_slices)
        res = root(f_min_spline,**kw)
        # get the parameter difference between each..
        u_x = res.x[0]
        # if the difference minimizing x and y is small enough, and
        # the minimization suceeded, and the minimization was small enough,
        # then we count this as a crossover
        succeeded = res.success
        small_value = (np.abs(res.fun) < tolerance).all()
        in_range = (u_x <= n_spline_slices)
        print(i,succeeded,small_value,in_range)
        if (succeeded and small_value and in_range):
            parameter_values_where_crossover.append(u_x)
    plt.close()
    for p in parameter_values_where_crossover:
        plt.plot(*PolymerTracing.evaluate_spline(p, tck=tck),color='r',marker='o')
    for s in spline_slices:
        plt.plot(*s)
    debug_idx = [21,23]
    for i_tmp in debug_idx:
        plt.plot(*spline_slices[i_tmp],color='r',linewidth=3)
        plt.plot(*spline_slices[min_dist_idx[i_tmp]],color='r',linewidth=3)

    plt.plot(x,y,'b.-')
    plt.show()


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
