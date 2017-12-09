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

sys.path.append("../../../../../../../../../../")
sys.path.append("../../")
from Util import IoUtil
from Research.Perkins.AnalysisUtil.Images import PolymerTracing
from GeneralUtil.python import CheckpointUtilities


def run(in_dir):
    """
    Args:
        in_dir: the input directory to operate on.  
    """
    input_dir =  IoUtil.data_dir(in_dir)
    cache_dir = IoUtil._traces_dir(in_dir)
    output_dir = IoUtil._ensemble_dir(in_dir)
    # just read in from the cache...
    objs_all = IoUtil.read_images(input_dir,cache_dir=cache_dir,force=False)
    kw = dict(min_m = 0,max_m = 125e-9)
    subset_funcs = [lambda o_tmp: o_tmp.dna_subset, 
                    lambda o_tmp: o_tmp.dna_protein_subset]
    data = []                        
    for subset_func in subset_funcs:
        polymer_info_obj = \
            PolymerTracing.ensemble_polymer_info(objs_all,
                                                 subset_func=subset_func,**kw)    
        data.append(polymer_info_obj)
    to_save = IoUtil.EnsembleObject(*data)
    CheckpointUtilities.lazy_save(output_dir + "ensemble.pkl",to_save)
    
if __name__ == "__main__":
    run(IoUtil.get_directory_command_line())
