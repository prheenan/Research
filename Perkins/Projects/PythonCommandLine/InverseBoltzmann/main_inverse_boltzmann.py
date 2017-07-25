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
    InverseBoltzmannUtil

def load_separation_wave(file_path):
    waves = PxpLoader.LoadAllWavesFromPxp(file_path)
    assert len(waves) == 1 , "Need exactly one wave"
    # POST: exactly one wave
    m_wave = waves[0]
    assert m_wave.Name().lower().endswidth("sep") , "Wave must end in 'sep'"
    # POST: this is a separation wave. return it 
    return m_wave


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
    interpolate_kwargs = dict(interpolation_factor=interpolation_factor)
    run_and_save_data(gaussian_stdev,extension,bins,interpolate_kwargs,out_file)
        

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
