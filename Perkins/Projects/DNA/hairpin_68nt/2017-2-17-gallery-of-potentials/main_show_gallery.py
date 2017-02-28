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

from scipy.interpolate import interp1d
from scipy.optimize import brute

def coth(x):
    return 1/np.tanh(x)

def fjc_extension(F,L0,Lp,kbT,K0):
    """

    See: 
    Wang, M. D., Yin, H., Landick, R., Gelles, J. & Block, S. M. 
    Stretching DNA with optical tweezers. Biophys J 72, 1335-1346 (1997)
    """
    return L0 * (coth( (2*F*Lp)/(kbT)) - kbT/(2*F*Lp)) * (1+F/K0)
    
def fjc_model_force(F,*args,**kwargs):
    force_range = [min(F),max(F)]
    return np.linspace(*force_range,num=F.size)
    
def fjc_model_ext_and_force(F,*args,**kwargs):
    force_modelled = fjc_model_force(F,*args,**kwargs)
    ext_modelled = fjc_extension(force_modelled,*args,**kwargs)
    return ext_modelled,force_modelled
    
def fjc_predicted_force_at_ext(separation,force,*args,**kwargs):
    ext_modelled,force_modelled = fjc_model_ext_and_force(force,*args,**kwargs)
    # interpolate back onto the grid we care about 
    interp = interp1d(ext_modelled,force_modelled,bounds_error=False,
                      fill_value="extrapolate")
    predicted_force = interp(separation)
    return predicted_force
    
def objective_l2(func_predict,true_values,*args):
    predicted_values = func_predict(*args)
    values = np.abs(predicted_values-true_values)**2
    to_ret =  sum(values)/sum(true_values**2)
    return to_ret

def brute_optimize(func_to_call,true_values,brute_dict=dict()):
    objective = lambda *args: objective_l2(func_to_call,true_values,*args)
    return brute(objective,disp=False,**brute_dict)
    
def fit_fjc_contour(separation,force,brute_dict=dict(),**kwargs):
    func = lambda *args: \
        fjc_predicted_force_at_ext(separation,force,*args,**kwargs)
    ret = brute_optimize(func,force,brute_dict=brute_dict)
    full_output_str = "full_output"
    if ( (full_output_str in brute_dict) and (brute_dict[full_output_str])):
        x0,fval,grid,Jout = ret
    else:
        x0 = ret
    model_x, model_y = fjc_model_ext_and_force(force,*x0,**kwargs)
    min_sep = min(separation)
    idx = np.where(model_x >= min_sep)
    return ret,model_x[idx],model_y[idx]
    
    
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
    force = False
    files,force_extension_curves = \
        FEC_Util.read_and_cache_pxp(directory,force=force)
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
    # get just the retracts
    curves_retract = [ fec.retract for fec in curves_split]
    curves_approach = [ fec.approach for fec in curves_split]
    # XXX debugging...
    ex = curves_retract[2]
    surface_idx = 1300
    max_idx = np.argmax(ex.Force)
    zeroed_force = ex.Force[surface_idx:max_idx]
    zeroed_separation = ex.Separation[surface_idx:max_idx]
    # ranges are for L0...
    max_sep = max(zeroed_separation)
    range_L0 = slice(0,max_sep,max_sep/50)
    brute_dict = dict(ranges=[range_L0])
    fjc_kwargs = dict(Lp=2.22e-9,kbT=4.1e-21,K0=1200e-12)
    ret,model_x,model_y = fit_fjc_contour(zeroed_separation,zeroed_force,
                                          brute_dict=brute_dict,**fjc_kwargs)
    print(ret)
    plt.plot(zeroed_separation,zeroed_force)
    plt.plot(model_x,model_y)
    plt.show()
    exit(1)
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
        save_name = "{:s}fec_{:d}_{:s}_{:s}.png".\
            format(out_dir,i,src_file,appr.Meta.Name)
        PlotUtilities.savefig(fig,save_name)
        
    
if __name__ == "__main__":
    run()
