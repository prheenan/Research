# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,ntpath
sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util import Analysis
from GeneralUtil.python import PlotUtilities


def fjc_extension(L0,F,Lp,kbT,K0):
    """

    See: 
    Wang, M. D., Yin, H., Landick, R., Gelles, J. & Block, S. M. 
    Stretching DNA with optical tweezers. Biophys J 72, 1335–1346 (1997)
    """
    return L0 * np.coth( (2*F*Lp)/(kbT) - kbT/(2*F*Lp)) * (1+F/K0)

def objective_l2(func_predict,true_values,*args):
    predicted_values = func_predict(*args)
    return np.abs(predicted_values-true_values)**2

def brute_optimize(func_to_call,brute_dict=None,*args,**kwargs):


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base_data_dir = FEC_Util.default_data_root()
    directory = base_data_dir + "4Patrick/CuratedData/DNA/Hairpin-68nt/" + \
       "2017-2-17-velocity-assay-single-attachments-no-adhesion" +  \
       "-1000nm-s-needs-to-be-redone/"
    out_dir = "./out/"
    files,force_extension_curves = \
        FEC_Util.read_and_cache_pxp(directory,force=False)
    force_extension_curves = [f for f in force_extension_curves]
    # split everything 
    curves_split= \
        [Analysis.split_FEC_by_meta(f) for f in force_extension_curves]
    # zero everything (of curves split, by reference)
    n_filter= lambda f: int(f.retract.Force.size*0.01)
    kwargs = dict(zero_separation_at_zero_force=True)
    [FEC_Util.zero_split_fec_approach_and_retract(f,NFilterPoints=n_filter(f),
                                                  **kwargs) 
     for f in curves_split]
    # XXX really awful  filtering
    curves_split = [c for i,c in enumerate(curves_split) 
                    if i not in [19,23,46,56,58,60,61,65,72]]
    # get just the retracts
    curves_retract = [ fec.retract for fec in curves_split]
    curves_approach = [ fec.approach for fec in curves_split]
    # plot as a 2-d histogram
    # filter to 1% of the size of the signal (should be loading-rate proof)
    histogram_filter_func = lambda o: int(np.ceil(o.Force.size * 0.01))
    fig = PlotUtilities.figure()
    max_separation_nanometers = 60
    FEC_Plot.heat_map_fec(curves_retract,n_filter_func= None,num_bins=(200,200),
                          separation_max=max_separation_nanometers)
    # XXX need to fix zeroing
    # limits in <nm,pN> for <x,y>
    plt.xlim([-10,60])
    plt.ylim([-30,60])
    n_base_pairs = 68
    # see (XXX need zotero citation)
    # http://www.sciencedirect.com/science/article/pii/S0378437112008771
    rise_nm_per_bp = 0.676
    contour_length_nm = n_base_pairs *rise_nm_per_bp
    style_contour = dict(color='w',linestyle='--',linewidth=4)
    label_contour = r"L$_0$" +"={:.1f} nm ({:d}bp)".format(contour_length_nm,
                                                           n_base_pairs)
    plt.axvline(contour_length_nm,label=label_contour,**style_contour)
    PlotUtilities.legend(loc='lower center',frameon=True,facecolor='w')
    PlotUtilities.savefig(fig,out_dir + "./out_histogram.png")   
    # plot each curve individually
    for i,(appr,retract) in enumerate(zip(curves_approach,curves_retract)):
        fig = PlotUtilities.figure()
        FEC_Plot.FEC_AlreadySplit(appr,retract,
                                  NFilterPoints=histogram_filter_func(retract))
        # use sensible units, <x,y> are <nm,pN>
        plt.xlim([-25,75])
        plt.ylim([-30,100])
        _,src_file = ntpath.split(appr.Meta.SourceFile)
        save_name = "{:s}fec_{:d}_{:s}.png".\
            format(out_dir,i,src_file,appr.Meta.Name)
        PlotUtilities.savefig(fig,save_name)
        
    
if __name__ == "__main__":
    run()
