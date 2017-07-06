# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,cProfile,os

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
from Research.Personal.EventDetection.Util import Analysis
from scipy.stats import norm
from FitUtil.EnergyLandscapes.Inverse_Boltzmann.Python.Code import \
    InverseBoltzmann

from scipy.interpolate import interp1d,griddata

class deconvolution_info:
    def __init__(self,ext_bins,ext_interp,interp_prob,P_q,P_q_interp,
                 p_k_interp):
        self.ext_bins = ext_bins
        self.ext_interp = ext_interp
        self.interp_prob = interp_prob
        self.P_q = P_q
        self.P_q_interp = P_q_interp
        self.p_k_interp = p_k_interp
        self.p_k = griddata(points=ext_interp,values=p_k_interp,xi=ext_bins)
        self.p_k /= np.trapz(y=self.p_k,x=self.ext_bins)
        self.free_energy_kT = -np.log(self.p_k)

def spring_const_plot(slices_safe,mean_safe_ext,std_safe_ext,pred):
    for i,s in enumerate(slices_safe):
        sep = s.Separation
        mean,std = mean_safe_ext[i],std_safe_ext[i]
        n,_,_ = plt.hist(sep,normed=True)
        plt.axvline(mean-std)
        plt.axvline(mean+std)
    plt.plot(mean_safe_ext,pred,color='r',linewidth=2)      
    PlotUtilities.lazyLabel("Extension (nm)","Standard Deviation (nm), PSF","")    

def probability_plot(inf):
    plt.subplot(2,1,1)
    plt.plot(inf.ext_bins,inf.p_k)
    plt.plot(inf.ext_bins,inf.P_q,'rp')
    plt.plot(inf.ext_interp,inf.p_k_interp,linestyle='--')
    plt.plot(inf.ext_interp,inf.P_q_interp,linestyle='--')
    PlotUtilities.lazyLabel("","PDF","")    
    plt.subplot(2,1,2)
    plt.plot(inf.ext_bins,inf.free_energy_kT)
    PlotUtilities.lazyLabel(r"Extension ($\AA$)","Free energy","")    

def deconvolution_plot(retract,slice_eq,slices,slices_safe,inf,
                       NFilterPoints,bins):
    plt.subplot(4,1,1)
    FEC_Plot._fec_base_plot(retract.Time,retract.Separation)
    plt.plot(slice_eq.Time,slice_eq.Separation,color='r',alpha=0.3)
    for s in slices:
        plt.axvline(s.Time[0]) 
        plt.axvline(s.Time[-1],linestyle='--')   
    plt.subplot(4,1,2)
    FEC_Plot._fec_base_plot(retract.Time,retract.Force,
                            n_filter_points=NFilterPoints)
    for s in slices:
        plt.axvline(s.Time[0])
    plt.plot(slice_eq.Time,slice_eq.Force,color='r',alpha=0.3)
    for s_safe in slices_safe:
        plt.plot(s_safe.Time,s_safe.Force,
                 color='g',alpha=0.3,linestyle='--')
    plt.subplot(4,1,3)
    plt.hist(slice_eq.Separation,normed=False,bins=bins)
    plt.subplot(4,1,4)
    plt.plot(inf.ext_bins*1e-10,inf.p_k)
    plt.plot(inf.ext_interp*1e-10,inf.P_q_interp,'r--')  
    plt.plot(inf.ext_interp * 1e-10,inf.p_k_interp,'b-')
    
def deconvolution(slice_eq,coeffs,bins):
    ext_eq = slice_eq.Separation * 1e10
    P_q,ext_bins = np.histogram(ext_eq,bins=bins,normed=True)
    # get rid of rightmost bin 
    ext_bins = ext_bins[:-1]
    # fit a spline to the bins to get a higher resolution histogram 
    interp_prob = interp1d(ext_bins,P_q,kind='cubic')
    ext_interp = np.linspace(min(ext_bins),max(ext_bins),endpoint=True,
                             num=bins*4)
    P_q_interp = interp_prob(ext_interp)
    P_q_interp /= np.trapz(y=P_q_interp,x=ext_interp)
    mean_ext_eq = np.mean(slice_eq.Separation)
    kwargs_deconv = dict(gaussian_stdev=np.polyval(coeffs,mean_ext_eq)*1e10,
                         extension=ext_interp,
                         P_q=P_q_interp)
    p_k_interp = InverseBoltzmann.gaussian_deconvolve(**kwargs_deconv)
    return deconvolution_info(ext_bins,ext_interp,interp_prob,P_q,P_q_interp,
                              p_k_interp)
                              
def get_slices(retract,NFilterPoints):                              
    fudge = 0.1
    delta_1 = 0.5
    offset_1 = 2.25
    n_1 = 4
    fudge = 0.01
    delta_2 = 2
    offset_2 = 4.259
    n_2 = 20
    fudge_2 = 0.1    
    delta_3 = 0.5
    offset_3 = 42.3
    time_slices = [ [offset_1 + delta_1 * i,offset_1 + delta_1 * (i+1)-fudge]
                    for i in range(n_1)] + \
                  [ [offset_2 + delta_2 * i,offset_2 + delta_2 * (i+1)-fudge_2]
                    for i in range(n_2)] + \
                  [ [offset_3 + delta_3 * i,offset_3 + delta_3 * (i+1)-fudge]
                    for i in range(4)]
    slices = [FEC_Util.slice_by_time(retract,*t)
               for t in time_slices]  
    # filter each slice independently 
    slices = [FEC_Util.GetFilteredForce(r,NFilterPoints=NFilterPoints)
              for r in slices]
    slice_idx_safe = [0,1,2,3,-4,-3,-2,-1]
    slices_safe = [slices[i] for i in slice_idx_safe]  
    return slices,slices_safe 
    
def analyze(example,out_dir):
    split_fec = Analysis.zero_and_split_force_extension_curve(example)
    retract = split_fec.retract
    dt = retract.Time[1] - retract.Time[0]
    NFilterPoints = int(np.ceil(5e-3/dt))
    slices,slices_safe  = get_slices(retract,NFilterPoints)
    mean_safe_ext = [np.mean(s.Separation) for s in slices_safe]
    std_safe_ext = [np.std(s.Separation) for s in slices_safe]
    coeffs,cov = np.polyfit(x=mean_safe_ext,y=std_safe_ext,full=False,cov=True,
                            deg=1)
    errors = np.sqrt(np.diag(cov))
    pred = np.polyval(coeffs,mean_safe_ext)
    pred_stderr = np.polyval(coeffs + errors,mean_safe_ext)
    # make the spring constant plot...
    fig = PlotUtilities.figure()
    spring_const_plot(slices_safe,mean_safe_ext,std_safe_ext,pred)
    PlotUtilities.savefig(fig,"{:s}spring.png".format(out_dir))
    # get a specific one for the equilibrium measurements 
    bins = 30
    for idx_eq in range(16,20):
        slice_eq = slices[idx_eq]
        inf = deconvolution(slice_eq,coeffs,bins=bins)
        # reinterpolate p_k_interp back onto the original grid 
        prob_savename = "{:s}probability{:d}.png".format(out_dir,idx_eq)
        fig = PlotUtilities.figure(figsize=(4,8))
        probability_plot(inf)
        PlotUtilities.savefig(fig,prob_savename)
        # make the deconvolution plot
        deconv_name = "{:s}deconvolution{:d}.png".format(out_dir,idx_eq)
        fig = PlotUtilities.figure(figsize=(4,8))
        deconvolution_plot(retract,slice_eq,slices,slices_safe,inf,
                           NFilterPoints,bins)
        PlotUtilities.savefig(fig,deconv_name)
        
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    abs_dir = "./"
    out_dir = "./out"
    cache_dir = "./cache/"
    examples = FEC_Util.\
        cache_individual_waves_in_directory(pxp_dir=abs_dir,force=False,
                                            cache_dir=cache_dir,limit=20)
    example = examples[-1]
    for example in examples:
        out_dir_tmp = "{:s}_{:s}/".format(out_dir,example.Meta.Name)
        GenUtilities.ensureDirExists(out_dir_tmp)
        analyze(example,out_dir_tmp)
        

    
if __name__ == "__main__":
    run()
