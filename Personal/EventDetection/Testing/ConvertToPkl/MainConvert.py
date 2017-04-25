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
sys.path.append("../../../../../")
from GeneralUtil.python import GenUtilities,CheckpointUtilities
from Research.Personal.EventDetection.Util import InputOutput



def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    directory_to_convert = InputOutput._base_curated() + \
                           "Positive/4nug2-devin/csv/"
    cache_directory = "./out/"
    GenUtilities.ensureDirExists(cache_directory)
    files = GenUtilities.getAllFiles(directory_to_convert,ext=".csv")
    for f in files:
        InputOutput.read_and_cache_file(f,cache_directory,
                                        has_events=True,force=True)

if __name__ == "__main__":
    run()
