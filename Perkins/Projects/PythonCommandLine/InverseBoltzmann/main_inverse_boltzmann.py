# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import sys
import matplotlib.pyplot as plt

import os, sys,traceback
sys.path.append('../../../../../')
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import GenUtilities
from IgorUtil.PythonAdapter import PxpLoader
import argparse
from FitUtil.EnergyLandscapes.Inverse_Boltzmann.Python.Code import \
    InverseBoltzmann

def load_separation_wave(file_path):
    waves = PxpLoader.LoadAllWavesFromPxp(file_path)
    assert len(waves) == 1 , "Need exactly one wave"
    # POST: exactly one wave
    m_wave = waves[0]
    assert m_wave.Name().lower().endswidth("sep") , "Wave must end in 'sep'"
    # POST: this is a separation wave. return it 
    return m_wave

def run_deconvolution(gaussian_stdev,extension,bins,
                      interpolate_kwargs = dict(),
                      deconvolve_common_kwargs=dict(p_0=None,
                                                    n_iters=300,
                                                    delta_tol=1e-9,
                                                    return_full=False,
                                                    r_0=1)):
    # get the extension distribution in whatever units the user gives us
    bins,P_q = \
        InverseBoltzmann.get_extension_bins_and_distribution(extension,
                                                             bins=bins)
    # XXX assume we know initial guess...
    p_0 = np.ones(P_q.size)
    sum_initial = sum(p_0)
    # get the normalized p_0
    p_0_normalized = p_0/np.trapz(y=p_0,x=bins)
    # determine what p_0 will then sum to
    p_0_sum = sum(p_0_normalized)
    # choose bins such that the sum is 1
    extension_factor = p_0_sum
    # XXX assume p0...
    extension_unitless = extension*extension_factor
    bins *= extension_factor
    gaussian_stdev *= extension_factor
    P_q /= np.trapz(y=P_q,x=bins)
    deconvolve_kwargs = dict(gaussian_stdev=gaussian_stdev,
                             extension_bins = bins,
                             P_q = P_q,
                             interpolate_kwargs=interpolate_kwargs,
                             **deconvolve_common_kwargs)
    interp_ext,interp_prob,deconv_interpolated_probability = \
        InverseBoltzmann.\
        interpolate_and_deconvolve_gaussian_psf(**deconvolve_kwargs)
    # convert the extensions back to their unit-less format, and renormalize the
    # probabilities so that they match up
    interp_ext = interp_ext * 1/extension_factor
    # 'raw' probability
    interp_prob /= np.trapz(x=interp_ext,y=interp_prob)
    # deconvolved probability 
    factor_deconv = np.trapz(x=interp_ext,y=deconv_interpolated_probability)
    deconv_interpolated_probability /= factor_deconv
    return interp_ext,interp_prob,deconv_interpolated_probability


def parse_and_run():
    parser = argparse.ArgumentParser(description='Inverse boltzmann of a .pxp ')
    common = dict(required=True)
    parser.add_argument('-number_of_bins', metavar='number_of_bins', 
                        type=int,help='number of approach/retract pairs',
                        **common)
    parser.add_argument('-interpolation_factor', 
                        metavar='interpolation_factor',
                        type=float,
                        help='force at which half the pop is folded/unfolded')
    help_gauss = "standard deviaiton of the (assumed gaussian) psf"
    parser.add_argument('-gaussian_stdev',metavar='gaussian_stdev',
                        type=float,help=help_gauss)
    parser.add_argument('-file_input',metavar="file_input",type=str,
                        help="path to the '.pxp' with the separation wave",
                        **common)
    parser.add_argument('-file_output',metavar="file_output",type=str,
                        help="path to output the associated data",**common)
    args = parser.parse_args()
    out_file = os.path.normpath(args.file_output)
    in_file = os.path.normpath(args.file_input)
    bins = args.number_of_bins
    gaussian_stdev = args.gaussian_stdev
    interpolation_factor = args.interpolation_factor
    data = load_separation_wave(in_file)
    # get the extension in ~ unitless for (it will return to 'normal' after)
    extension = data.dataY
    interpolate_kwargs=dict(interpolation_factor=interpolation_factor)
    args = run_deconvolution(gaussian_stdev,extension,bins,
                             interpolate_kwargs=interpolate_kwargs)
    interp_ext,interp_prob,deconv_interpolated_probability = args
        

def run():
    # change to this scripts path
    path = os.path.abspath(os.path.dirname(__file__))
    os.chdir(path)
    try:
        parse_and_run()
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lines = traceback.format_exception(exc_type, exc_value, 
                                           exc_traceback)
        # Log it or whatever here
        str_out =''.join('!! ' + line for line in lines)
        print(str_out)
        
if __name__ == "__main__":
    run()
