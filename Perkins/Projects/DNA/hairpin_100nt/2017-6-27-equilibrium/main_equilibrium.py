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
from GeneralUtil.python import PlotUtilities,CheckpointUtilities,GenUtilities
from Research.Personal.EventDetection.Util import Analysis
from scipy.stats import norm
from FitUtil.EnergyLandscapes.Inverse_Boltzmann.Python.Code import \
    InverseBoltzmann,InverseBoltzmannUtil

from scipy.interpolate import interp1d,griddata

class deconvolution_info:
    def __init__(self,ext_bins,ext_interp,interp_prob,P_q,
                 p_k_interp):
        self.ext_bins = ext_bins
        self.ext_interp = ext_interp
        self.interp_prob = interp_prob
        self.P_q = P_q
        self.p_k_interp = p_k_interp
        self.p_k = griddata(points=ext_interp,values=p_k_interp,xi=ext_bins)
        self.p_k /= np.trapz(y=self.p_k,x=self.ext_bins)
        self.free_energy_kT = -np.log(self.p_k)
        self.convolved_free_energy_kT = -np.log(self.P_q)

def spring_const_plot(slices_safe,mean_safe_ext,std_safe_ext,pred):
    plt.plot(mean_safe_ext,std_safe_ext,'ro-',linewidth=2)      
    plt.plot(mean_safe_ext,pred,'b--')
    PlotUtilities.lazyLabel("Extension (nm)","Standard Deviation (nm), PSF","")    
def histogram_plot(s):
    sep = s.Separation * 1e9
    force_to_plot = s.Force*1e12    
    plt.subplot(2,2,1)
    plt.plot(s.Time,sep,color='g',linewidth=0.2)
    PlotUtilities.lazyLabel("","Separation (nm)","")        
    plt.subplot(2,2,2)
    plt.hist(sep,orientation='horizontal',color='g')
    mean,stdev = np.mean(sep),np.std(sep)
    plt.axhline(mean-stdev)
    plt.axhline(mean+stdev,label=r"$\sigma$=" + "{:.2g}nm".format(stdev))
    PlotUtilities.lazyLabel("","Count","")            
    plt.subplot(2,2,3)
    plt.plot(s.Time,force_to_plot,color='r',linewidth=0.25)
    PlotUtilities.lazyLabel("Time (S)","Force (pN)","")            
    plt.subplot(2,2,4)
    plt.hist(force_to_plot,orientation='horizontal',color='r')
    mean,stdev = np.mean(force_to_plot),np.std(force_to_plot)
    plt.axhline(mean-stdev)
    plt.axhline(mean+stdev,label=r"$\sigma$=" + "{:.2g}pN".format(stdev))
    PlotUtilities.lazyLabel("Count","","")          

def probability_plot(inf,slice_eq):
    plt.subplot(3,1,1)
    plt.plot(inf.ext_bins,inf.p_k)
    plt.plot(inf.ext_bins,inf.P_q,'r-')
    plt.plot(inf.ext_interp,inf.p_k_interp,linestyle='--')
    PlotUtilities.lazyLabel("","PDF","")    
    plt.subplot(3,1,2)
    plt.plot(inf.ext_bins,inf.free_energy_kT)
    plt.plot(inf.ext_bins,inf.convolved_free_energy_kT)
    PlotUtilities.lazyLabel(r"Extension ($\AA$)","Energy at <F>","")        
    plt.subplot(3,1,3)
    bins_relative = inf.ext_bins + min(inf.ext_bins)
    tilt_energy  = (np.mean(slice_eq.Force) * bins_relative)/4.1e-21
    plt.plot(inf.ext_bins,inf.free_energy_kT-tilt_energy)
    PlotUtilities.lazyLabel(r"Extension ($\AA$)","Free energy","")    

def deconvolution_plot(retract,slice_eq,slices,inf,
                       NFilterPoints,bins):
    sep_limits = [min(inf.ext_bins),max(inf.ext_bins)]                      
    time_limits = [min(slice_eq.Time),max(slice_eq.Time)]
    retract_sliced_tmp = FEC_Util.slice_by_time(retract,*time_limits)
    plt.subplot(4,1,1)
    plt.plot(slice_eq.Time,slice_eq.Separation*1e10,color='r',alpha=0.3)
    plt.plot(retract_sliced_tmp.Time,retract_sliced_tmp.Separation*1e10,color='k',
             alpha=0.3)
    plt.xlim(time_limits)
    PlotUtilities.lazyLabel(r"",r"Separation ($\AA$)","")    
    plt.subplot(4,1,2)
    plt.plot(slice_eq.Time,slice_eq.Force*1e12,color='r',alpha=0.3)
    plt.plot(retract_sliced_tmp.Time,retract_sliced_tmp.Force*1e12,color='k',
             alpha=0.3)    
    plt.xlim(time_limits)             
    PlotUtilities.lazyLabel(r"Time (s)","Force (pN)","")        
    plt.subplot(4,1,3)
    plt.hist(slice_eq.Separation,normed=False,bins=bins)
    plt.xlim(sep_limits)    
    PlotUtilities.lazyLabel("","Count","")            
    plt.subplot(4,1,4)
    plt.plot(inf.ext_bins,inf.p_k)
    plt.plot(inf.ext_interp,inf.p_k_interp,'b-')
    plt.xlim(sep_limits)        
    PlotUtilities.lazyLabel(r"Separation ($\AA$)","PDF","")            
    
    
def deconvolution(slice_eq,coeffs,bins):
    ext_eq = slice_eq.Separation
    # get the 'raw' distributions
    ext_bins,P_q = InverseBoltzmannUtil.\
        get_extension_bins_and_distribution(ext_eq,bins=bins)  
    mean_ext_eq = np.mean(ext_eq)            
    # XXX rough estimate for stdev
    gaussian_stdev=np.polyval(coeffs,mean_ext_eq)
    # do the actual deconvolution 
    args = InverseBoltzmannUtil.extension_deconvolution(gaussian_stdev,
                                                        ext_eq,bins)
    interp_ext,interp_prob,p_k_interp = args
    return deconvolution_info(ext_bins,interp_ext,interp_prob,P_q,
                              p_k_interp)
                              
def get_slices(retract,NFilterPoints):                              
    fudge = 0.1
    delta_1 = 0.5
    offset_1 = 2.25
    n_1 = 4
    fudge = 0.01
    delta_2 = 2
    offset_2 = 4.259
    n_2 = 19
    fudge_2 = 0.1    
    delta_3 = 0.5
    offset_3 = 42.27
    fudge_3 = 0.1
    time_slices = [ [offset_1 + (delta_1) * i + fudge,
                     offset_1 + (delta_1) * (i+1)]
                    for i in range(n_1)] + \
                  [ [offset_2 + (delta_2) * i + fudge,
                    offset_2 + delta_2 * (i+1)]
                    for i in range(n_2)] + \
                  [ [offset_3 + (delta_3)* i,
                    offset_3 + delta_3 * (i+1)-fudge]
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
    for i,s in enumerate(slices):
        fig = PlotUtilities.figure()
        histogram_plot(s)
        PlotUtilities.savefig(fig,"{:s}_hist_{:d}.png".format(out_dir,i))   
    # get a specific one for the equilibrium measurements 
    bins = 100
    for idx_eq in range(0,24):
        slice_eq = slices[idx_eq]
        inf = deconvolution(slice_eq,coeffs,bins=bins)
        # POST: deconvolution worked
        # reinterpolate p_k_interp back onto the original grid 
        prob_savename = "{:s}probability{:d}.png".format(out_dir,idx_eq)
        fig = PlotUtilities.figure(figsize=(4,8))
        probability_plot(inf,slice_eq)
        PlotUtilities.savefig(fig,prob_savename)
        # make the deconvolution plot
        deconv_name = "{:s}deconvolution{:d}.png".format(out_dir,idx_eq)
        fig = PlotUtilities.figure(figsize=(4,8))
        deconvolution_plot(retract,slice_eq,slices,inf,
                           NFilterPoints,bins)
        PlotUtilities.savefig(fig,deconv_name)
        CheckpointUtilities.lazy_save(prob_savename +".pkl",inf)
        
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    abs_dir = "./"
    out_dir = "./out/"
    cache_dir = "./cache/"
    examples = FEC_Util.\
        cache_individual_waves_in_directory(pxp_dir=abs_dir,force=False,
                                            cache_dir=cache_dir,limit=20)
    for example in examples:
        out_dir_tmp = "{:s}{:s}/".format(out_dir,example.Meta.Name)
        GenUtilities.ensureDirExists(out_dir_tmp)
        analyze(example,out_dir_tmp)
        

    
if __name__ == "__main__":
    run()
