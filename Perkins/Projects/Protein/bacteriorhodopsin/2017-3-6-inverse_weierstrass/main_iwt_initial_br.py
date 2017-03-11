# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

def hao_grouping_function(str_v):
    pattern = r"""
               (\D+) # type (eg: ext,force)
               (\d+)      # id (e.g. 1131)
               (\D+)        # anything else, who cares
               """
    match = re.match(pattern,str_v,re.VERBOSE)
    assert match is not None , "Whoops! Got a bad string: {:s}".format(str_v)
    ending,id,preamble = match.groups()
    # convert ext to sep
    ending = ending if ending != "ext" else "sep"
    return preamble,id,ending

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = FEC_Util.default_data_root()
    relative_data_base = "/HighSpeedPull/FEC_30nms/"
    relative_data_dir = "161011_FEC_raw/"
    absolute_data_dir = base + relative_data_base + relative_data_dir
    # XXX use Haos...
    absolute_data_dir = base + "4Patrick/Scratch/Tmp_Data_Scratch/"
    dict_kwargs = dict(cache_directory="./",
                       in_directory=absolute_data_dir,
                       grouping_function = hao_grouping_function,
                       limit=3,force=True)
    data = FEC_Util.cache_ibw_directory(**dict_kwargs)  

if __name__ == "__main__":
    run()
