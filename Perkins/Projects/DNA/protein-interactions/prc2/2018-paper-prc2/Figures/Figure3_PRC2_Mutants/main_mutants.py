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

sys.path.append("../../../../../../../../../")
sys.path.append("../")
from Util import IoUtil
from Research.Perkins.AnalysisUtil.Images import PolymerTracing
from GeneralUtil.python import CheckpointUtilities,PlotUtilities

hist_default_dict = dict(normed=True,linewidth=0,color='b',alpha=0.3)

def contour_length_plot(ensemble,hist_kwargs=hist_default_dict):
    L0_nm = np.array(ensemble.L0_m) * 1e9
    hist = plt.hist(L0_nm,label="$N$={:d}".format(L0_nm.size),**hist_kwargs)
    mean_v = np.mean(L0_nm)
    percentiles = np.percentile(L0_nm,[5,95])
    f_round = lambda x : int(np.round(x,-1))
    mean_nm_rounded = f_round(mean_v)
    percentile_label = "[5%,95%] = [{:d},{:d}] nm".\
        format(*[f_round(p) for p in percentiles])
    plt.axvline(mean_v,color='r',linestyle='--',
                label="$\mu$ = {:d} nm".format(mean_nm_rounded))
    for i,p in enumerate(percentiles):
        label_tmp = percentile_label if i == 0 else ""
        plt.axvline(p,label=label_tmp,color='g',linestyle=":")
    PlotUtilities.lazyLabel("$L_0$ (nm)","P (1/nm)","",frameon=True)
    return hist 
    
def persistence_length_plot(ensemble,fit_bounds_nm=[0,60]):
    minus_log_mean_angle = -np.log(ensemble.cos_angle_binned)
    max_length_nm = max(fit_bounds_nm)
    L_binned_plot = ensemble.L_binned * 1e9
    fit_idx = np.where( (L_binned_plot >= min(fit_bounds_nm)) & \
                        (L_binned_plot <= max_length_nm) )
    assert fit_idx[0].size > 0 , "Bounds incorrect "
    fit_x = L_binned_plot[fit_idx]
    fit_y = minus_log_mean_angle[fit_idx]
    coeffs,cov = np.polyfit(x=fit_x,y=fit_y,deg=1,cov=True)
    # square root in the diagonal is the uncertainty in the fitting parameters
    coeffs_error = np.sqrt(np.diag(cov))
    #     <Cos<Angle>>    = e(-L/(L_p))
    # -Log(<Cos<Angle>>)  = L/(L_p)
    # so on a plot of -Log(<Cos<Angle>>) vs L, the offset should be zero,
    # the slope should be 
    # slope = 1/(L_p)
    Lp = 1/coeffs[0]
    # L_p = 1/c0
    # sigma_(L_p) = (1/(c0^2)) * sigma_c0 
    #             = Lp^2 * sigma_c0
    L_p_error = (Lp**2) * coeffs_error[0]
    fit_label = "$L_p$ = {:.1f} $\pm$ {:.1f} nm".format(Lp,L_p_error)
    interp_x = np.linspace(min(fit_x),max(fit_x),endpoint=True,num=fit_x.size*5)
    interp_y = np.polyval(coeffs,x=interp_x)
    plt.plot(L_binned_plot,minus_log_mean_angle)
    plt.plot(interp_x,interp_y,'r--',label=fit_label)
    plt.axvspan(max_length_nm,max(plt.xlim()),color='k',alpha=0.3,
                linewidth=0,
                label="Not fit")
    PlotUtilities.lazyLabel("Length $L$ (nm)",r"-$\log{<\cos(\theta)>}$","")


def plot_specfic_directory(in_dir):
    """
    Args:
        in_dir: the input directory to operate on.  
    """
    input_dir =  IoUtil._ensemble_dir(in_dir)
    obj = CheckpointUtilities.lazy_multi_load(input_dir)
    output_dir = IoUtil._plot_dir(in_dir)
    xlim_L0_nm = [0,1e3]
    xlim_L_nm = [0,110]
    ylim_probability = [0,1.2e-2]
    ylim_minus_log_cos = [0,5]
    assert len(obj) > 0 , "No data saved in {:s}".format(input_dir)
    ensemble = obj[0]
    output_tuples = [ [ensemble.dna_only,"DNA Only"],
                      [ensemble.dna_plus_protein,"DNA + Protein"]]
    for data,name in output_tuples:
        # make a figure of the dna...
        fig = PlotUtilities.figure((7.5,2.5))
        ax_L0 = plt.subplot(1,2,1)
        plt.xlim(xlim_L0_nm)
        plt.ylim(ylim_probability)       
        contour_length_plot(data)        
        ax_Lp = plt.subplot(1,2,2)
        plt.xlim(xlim_L_nm)
        plt.ylim(ylim_minus_log_cos)
        persistence_length_plot(data)
        PlotUtilities.savefig(fig,output_dir + "polymer{:s}.png".format(name))

def run():
    directories = ["../../Data/601-widom/5mer/1x/",
                   "../../Data/designed/at-rich/1x/",
                   "../../Data/designed/gc-rich/1x/"]
    in_dir = [IoUtil._ensemble_dir(d) for d in directories]
    xlim_L0_nm = [0,1e3]
    ensembles = [CheckpointUtilities.lazy_multi_load(d) for d in in_dir]
    hist_arr = []
    fig = PlotUtilities.figure((3.5,7))
    ax_arr = []
    n = len(ensembles)
    hist_kw = dict(hist_default_dict)
    binsize_nm = 50
    max_L0 = max(xlim_L0_nm)
    N = max_L0/binsize_nm
    hist_kw['bins'] = np.linspace(0,max_L0,num=20,endpoint=True)
    colors = ['g','r','b','k']
    titles = ["12x601","GC-poor","GC-rich center"]
    for i,e in enumerate(ensembles):
        ax_tmp = plt.subplot(n,1,(i+1))
        hist_kw_tmp = dict(**hist_kw)
        hist_kw_tmp['color'] = colors[i]
        hist_tmp = contour_length_plot(e[0].dna_plus_protein,
                                       hist_kwargs=hist_kw_tmp)
        hist_arr.append(hist_tmp)
        ax_arr.append(ax_tmp)
        PlotUtilities.title(titles[i],color=colors[i])
        if (i != n-1):
            PlotUtilities.xlabel("")
            PlotUtilities.no_x_label(ax_tmp)
    max_n = [max(h[0]) for h in hist_arr]
    ylim_P = [0,max(max_n)]
    for a in ax_arr:
        a.set_xlim(xlim_L0_nm)
        a.set_ylim(ylim_P)
        PlotUtilities.tom_ticks(a,change_x=False,num_major=3)
    PlotUtilities.savefig(fig,"./out.png")
        

if __name__ == "__main__":
    run()

