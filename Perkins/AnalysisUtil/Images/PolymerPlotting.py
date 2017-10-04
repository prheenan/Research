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
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.Images import PolymerTracing

def plot_single_worm_object(tmp_fit):
    plt.subplot(2,1,1)
    plt.plot(*tmp_fit.fit_xy,color='g',marker='o',linewidth=0.5)
    plt.subplot(2,1,2)
    x,y = tmp_fit.x_y_rel
    plt.plot(x,y,color='r',linewidth=2)
    plt.imshow(tmp_fit.image_threshold)

def plot_angle_information(polymer_info_obj):
    L_binned = polymer_info_obj.L_binned
    log_mean_angle = -np.log(polymer_info_obj.cos_angle_binned)
    x,t1,t2,t3,t4 = PolymerTracing.theta_stats(polymer_info_obj,n_bins=50)
    log_t2 = np.log(t2)
    log_x = np.log(x)
    idx_min = np.argmin(abs(log_t2))
    log_x_min = log_x[idx_min]
    Lp_m = np.exp(log_x_min)
    to_x = lambda x: x*1e9
    plt.subplot(2,1,1)
    plt.plot(log_x,log_t2)
    plt.axvline(log_x_min,label=r"$L_p=$" + "{:.2g} nm".format(Lp_m*1e9))
    PlotUtilities.lazyLabel("Log(x) (log(nm))",
                            r"$\log(\theta^2)$","")
    plt.subplot(2,1,2)
    plt.plot(to_x(L_binned),log_mean_angle)
    PlotUtilities.lazyLabel("Extension",r"$\log{<\cos(\theta)>}$","")

