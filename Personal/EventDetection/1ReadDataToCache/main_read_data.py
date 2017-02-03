# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os

sys.path.append("../../../../")
from GeneralUtil.python import GenUtilities,CheckpointUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util


def read_and_cache_file(cache_directory,file_path,has_events=False,force=False):
    file_name = os.path.basename(file_path)
    cache_file = cache_directory + file_name+ ".pkl"
    func_to_call = FEC_Util.read_time_sep_force_from_csv
    return CheckpointUtilities.getCheckpoint(cache_file,func_to_call,force,
                                             file_path,has_events=has_events)

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base_directory = "/Volumes/group/4Patrick/CuratedData/Masters_CSCI/"
    positives_directory = base_directory + "Positive/"
    cache_directory = base_directory + "cache/"
    # relative directories are given by <directory,data type>
    # where data type is something like dna-650nm...
    relative_directories = \
            [["CircularDNA/ovh2.0-1p9kbp/ForceExtensionCurves/RawPxpFiles",
              "dna-650nm"]]
    absolute_positives = [positives_directory + r[0] 
                          for r in relative_directories]
    force = True
    # get the positive events
    for i,dir_v in enumerate(absolute_positives):
        # get the label for this dataset.
        label = relative_directories[i][1]
        all_files = GenUtilities.getAllFiles(dir_v,ext=".csv")
        for f in all_files:
            # XXX save this by label...
            read_and_cache_file(cache_directory,f,has_events=True,force=force)
    # get the negative events
    # XXX 


if __name__ == "__main__":
    run()
