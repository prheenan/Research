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

def plot_angle_information(polymer_info_obj,max_length_nm=75):
    L_binned = polymer_info_obj.L_binned
    minus_log_mean_angle = -np.log(polymer_info_obj.cos_angle_binned)
    # get the stdev in the bins 
    log_cos_sq = np.log(polymer_info_obj.cos_angle)**2
    L_binned_sq,log_angle_sq_binned = \
        PolymerTracing._binned_stat(x=polymer_info_obj.L_m,
                                    y=log_cos_sq,
                                    n_bins=L_binned.size)
    np.testing.assert_allclose(L_binned_sq,L_binned)                                    
    variance_log_cos = log_angle_sq_binned - minus_log_mean_angle**2
    std_log_cos = np.sqrt(variance_log_cos)
    x,t1,t2,t3,t4 = PolymerTracing.theta_stats(polymer_info_obj,n_bins=50)
    log_t2 = np.log(t2)
    log_x = np.log(x)
    idx_min = np.argmin(abs(log_t2))
    log_x_min = log_x[idx_min]
    Lp_m = np.exp(log_x_min)
    to_x = lambda x: x*1e9
    L_binned_plot = to_x(L_binned)
    fit_idx = np.where(L_binned_plot <= max_length_nm)[0]
    assert fit_idx.size > 0 , "Need positive lengths..."
    fit_x = L_binned_plot[fit_idx]
    fit_y = minus_log_mean_angle[fit_idx]
    interp_x = np.linspace(L_binned_plot[0],L_binned_plot[-1],endpoint=True,
                           num=int(L_binned_plot.size*10))
    coeffs,cov = np.polyfit(x=fit_x,y=fit_y,deg=1,cov=True)
    # square root in the diagonal is the uncertainty in the fitting parameters
    coeffs_error = np.sqrt(np.diag(cov))
    Lp = 1/coeffs[0]
    # L_p = 1/c0
    # sigma_(L_p) = (1/(c0^2)) * sigma_c0 
    #             = Lp^2 * sigma_c0
    L_p_error = (Lp**2) * coeffs_error[0]
    fit_label = "$L_p$ = {:.1f} $\pm$ {:.1f} nm".format(Lp,L_p_error)
    interp_y = np.polyval(coeffs,x=interp_x)
    plt.subplot(1,1,1)
    plt.plot(L_binned_plot,minus_log_mean_angle)
    plt.plot(interp_x,interp_y,'r--',label=fit_label)
    upper = minus_log_mean_angle + std_log_cos
    lower = minus_log_mean_angle - std_log_cos
    plt.fill_between(x=L_binned_plot,y1=lower,y2=upper,alpha=0.3,color='g')
    plt.axvspan(max_length_nm,max(plt.xlim()),color='k',alpha=0.3,
                linewidth=0,
                label="Not fit")
    PlotUtilities.lazyLabel("Length $L$ (nm)",r"$\log{<\cos(\theta)>}$","")

