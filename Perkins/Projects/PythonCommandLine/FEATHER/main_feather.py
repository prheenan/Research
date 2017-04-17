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
import argparse
from Research.Personal.EventDetection._2SplineEventDetector import Detector

def write_and_close(string):
    raise RuntimeError(string)

def parse_and_run():
    description = 'Predict event locations in a .pxp '
    parser = argparse.ArgumentParser(description=description)
    common = dict(required=True)
    parser.add_argument('-threshold', metavar='threshold', 
                        type=float,help='probability threshold (between 0 and 1)',
                        **common)
    parser.add_argument('-file_input',metavar="file_input",type=str,
                        help="path to the '.pxp' with the force, separation",
                        **common)
    parser.add_argument('-file_output',metavar="file_output",type=str,
                        help="path to output the associated data",**common)
    args = parser.parse_args()
    out_file = os.path.normpath(args.file_output)
    in_file = os.path.normpath(args.file_input)
    threshold = args.threshold
    if (not GenUtilities.isfile(in_file)):
        write_and_close("File {:s} doesn't exist".format(in_file))
    # # POST: input file exists
    # go ahead and read it
    RawData = FEC_Util.ReadInData(in_file,Limit=1)
    # POST: file read sucessfully. should just have the one
    if (not len(RawData) == 1):
        write_and_close("Need exactly one Force/Separation".\
                        format(in_file))
    # POST: have just one. Go ahead and using FEATHER to predict the locations
    example = RawData[0]
    event_indices = Detector.predict(example,threshold=threshold)
    print(event_indices)
    # done with the log file...
    np.savetxt(fname=out_file,delimiter=",",newline="\n",
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
