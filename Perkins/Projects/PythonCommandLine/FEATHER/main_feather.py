# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import sys
import matplotlib.pyplot as plt
import os, sys,traceback
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
import argparse,h5py
from Research.Personal.EventDetection._2SplineEventDetector import Detector

def write_and_close(string):
    raise RuntimeError(string)

def read_matlab_file_into_fec(input_file):
    """
    Reads a matlab file into a force extension curve

    Args:
        input_file: '.mat' file, formatted like -v7.3
    Returns:
        tuple of time,separation,force
    """
    f = h5py.File(input_file,'r') 
    get = lambda x: f[x].value
    # get the FEC data
    time = get('time')
    separation = get('separation')
    force = get('force')
    return time,separation,force
    

def get_force_extension_curve(in_file,**kwargs):
    """
    given an input file and meta information, returns the associated force 
    extension curve

    Args:
         input_file: .pxp or .m file 
    Returns:
         force extension curve object which FEATHER can use
    """
    if (not GenUtilities.isfile(in_file)):
        write_and_close("File {:s} doesn't exist".format(in_file))
    # # POST: input file exists
    # go ahead and read it
    if (in_file.endswith(".pxp")):
        RawData = FEC_Util.ReadInData(in_file,Limit=1)
        # POST: file read sucessfully. should just have the one
        if (not len(RawData) == 1):
            write_and_close("Need exactly one Force/Separation".\
                            format(in_file))
        # POST: have one. Go ahead and use FEATHER to predict the locations
        to_ret = RawData[0]
    elif (in_file.endswidth(".mat") or in_file.endswidth(".m")):
        time,separation,force = read_matlab_file_into_fec(input_file)
        data = TimeSepForceObj.data_obj_by_columns_and_dict(time=time,
                                                            sep=separation,
                                                            force=force,
                                                            meta_dict=meta_dict)
        to_ret = TimeSepForceObj.TimeSepForceObj()
        to_ret.LowResData = data
    return to_ret 



def parse_and_run():
    description = 'Predict event locations in a data file'
    parser = argparse.ArgumentParser(description=description)
    common = dict(required=True)
    # # feathers options
    parser.add_argument('-tau', metavar='tau', 
                        type=float,help='tau fraction of curve (0,1)',
                        required=False,default=1e-2)
    parser.add_argument('-threshold', metavar='threshold', 
                        type=float,help='probability threshold (0,1)',
                        **common)
    # # 'meta' variables
    parser.add_argument('-spring_constant', metavar='spring_constant', 
                        type=float,help='spring constant of the probe',
                        **common)
    parser.add_argument('-trigger_time', metavar='trigger_time', 
                        type=float,help='time at which approach ends',
                        **common)
    parser.add_argument('-dwell_time', metavar='dwell_time', 
                        type=float,
                        help='time between end of approach and retract start',
                        **common)
    # path to the file
    parser.add_argument('-file_input',metavar="file_input",type=str,
                        help="path to the force-extension curve file",
                        **common)
    parser.add_argument('-file_output',metavar="file_output",type=str,
                        help="path to output the associated data",**common)
    args = parser.parse_args()
    out_file = os.path.normpath(args.file_output)
    in_file = os.path.normpath(args.file_input)
    threshold = args.threshold
    tau = args.tau
    assert tau > 0 , "FEATHER Tau must be greater than 0"
    assert threshold > 0 , "FEATHER threshold must be greater than 0"
    example = get_force_extension_curve(in_file,
                                        K=args.spring_constant,
                                        DwellTime=args.dwell_time,
                                        TriggerTime=args.trigger_time,
                                        # set these to one; aren't interested
                                        # in volts (feather works with FECs)
                                        DwellSetting=1,
                                        Invols=1)
    event_indices = Detector.predict(example,threshold=threshold,
                                     add_offsets=True,tau_fraction=tau)
    # done with the log file...
    np.savetxt(fname=out_file,delimiter=",",newline="\n",fmt="%d",
               header="(C) PRH 2017\nEvent Indices",
               X=event_indices)

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
