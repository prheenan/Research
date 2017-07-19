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

sys.path.append("../../../../../../../")
from Research.Personal.EventDetection.Util import InputOutput

def run():
    """
    This file converts new .pxp files into csv files; useful if we have new edge
    cases coming in 
    """
    in_dir = "./data/"
    out_dir = "./out/"
    InputOutput.output_waves_in_directory_to_csv_files(in_dir,out_dir)

if __name__ == "__main__":
    run()
