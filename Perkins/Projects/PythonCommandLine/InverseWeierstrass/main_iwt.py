# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import sys

import os, sys,traceback
import matplotlib.pyplot as plt
# change to this scripts path
path = os.path.abspath(os.path.dirname(__file__))
os.chdir(path)
sys.path.append('../../../../../')
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import GenUtilities
from IgorUtil.PythonAdapter import PxpLoader
import argparse

def write_and_close(string):
    raise RuntimeError(string)


def parse_and_run():
    parser = argparse.ArgumentParser(description='IWT of a .pxp ')
    common = dict(required=True)
    parser.add_argument('-number_of_pairs', metavar='number_of_pairs', 
                        type=int,help='number of approach/retract pairs',
                        **common)
    parser.add_argument('-flip_forces',metavar="flip_forces",type=int,
                        help="if true, multiplies all the forces by -1",
                        **common)
    parser.add_argument('-number_of_bins', metavar='number_of_bins', 
                        type=int,help='number of separation bins',**common)
    parser.add_argument('-f_one_half', metavar='f_one_half', type=float,
                        help='force at which half the pop is folded/unfolded',
                        **common)
    help_vel = '[0,1] of the separation vs time to fit for the velocity'
    parser.add_argument('-fraction_velocity_fit', 
                        metavar='fraction_velocity_fit', type=float,
                        help=help_vel,**common)
    parser.add_argument('-file_input',metavar="file_input",type=str,
                        help="path to the '.pxp' with the force, separation",
                        **common)
    parser.add_argument('-file_output',metavar="file_output",type=str,
                        help="path to output the associated data",**common)
    vel_help = "optional manually-specified velocity (m/s). If this is" + \
               " present, then it is used to determing the velocity instead" +\
               " of fraction_velocity_fit"
    parser.add_argument('-velocity',metavar="velocity",type=float,default=0,
                        help=vel_help,required=False)
    args = parser.parse_args()
    out_file = os.path.normpath(args.file_output)
    in_file = os.path.normpath(args.file_input)
    flip_forces = args.flip_forces
    number_of_pairs = args.number_of_pairs
    f_one_half = args.f_one_half
    if (not GenUtilities.isfile(in_file)):
        write_and_close("File {:s} doesn't exist".format(in_file))
    # # POST: input file exists
    # go ahead and read it
    validation_function = PxpLoader.valid_fec_allow_endings
    RawData = IWT_Util.ReadInAllFiles([in_file],Limit=1,
                                      ValidFunc=validation_function)
    # POST: file read sucessfully. should just have the one
    if (not len(RawData) == 1):
        write_and_close("Need exactly one Force/Separation".\
                        format(in_file))
    # POST: have just one. Go ahead and break it up
    fraction_for_vel = args.fraction_velocity_fit    
    Data= RawData[0]
    unfold,refold = IWT_Util.\
        get_unfold_and_refold_objects(Data,
                                      number_of_pairs=number_of_pairs,
                                      flip_forces=flip_forces,
                                      fraction_for_vel=fraction_for_vel)
    # if the velocity was there, use it.
    velocity = args.velocity
    if ( (velocity is not None) and velocity > 0):
        for un,re in zip(unfold,refold):
            # keep the offsets, reset the velocites
            un.SetOffsetAndVelocity(un.Offset,velocity)
            re.SetOffsetAndVelocity(re.Offset,velocity * -1)
    # POST: have the unfolding and refolding objects, get the energy landscape
    num_bins = args.number_of_bins
    LandscapeObj =  InverseWeierstrass.\
        FreeEnergyAtZeroForce(unfold,NumBins=num_bins,RefoldingObjs=refold)
    # get the distance to the transition state etc
    all_landscape = [-np.inf,np.inf]    
    Obj =  IWT_Util.TiltedLandscape(LandscapeObj,f_one_half_N=f_one_half)
    # write out the file we need
    extension_meters = Obj.landscape_ext_nm/1e9
    landscape_joules = Obj.Landscape_kT * Obj.kT
    landscape_tilted_joules = Obj.Tilted_kT * Obj.kT
    data = np.array((extension_meters,landscape_joules,
                     landscape_tilted_joules))
    # done with the log file...
    np.savetxt(fname=out_file,delimiter=",",newline="\n",
               header="(C) PRH 2017\n extension(m),landscape(J),tilted(J)",
               X=data.T)


def run():
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
