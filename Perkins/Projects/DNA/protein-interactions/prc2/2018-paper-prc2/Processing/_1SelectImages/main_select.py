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

sys.path.append("../")
sys.path.append("../../../../../../../../../")
from Research.Perkins.AnalysisUtil.Images import ImageUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

from GeneralUtil.python import CheckpointUtilities
from Util import Processing

def filter(images,expected_length_m):
    cache_dir = Processing.cache_select_images()
    for i in images:
        if (abs(i.range_meters - 2e-6) < 1e-12):
            yield i

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base_dir = Processing.cache_raw_images()
    expected_length_m = 2e-6
    cache_dir = Processing.cache_select_images()
    images = CheckpointUtilities.lazy_multi_load(base_dir)
    load_func = lambda: filter(images,expected_length_m)
    e = CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_func,
                                       name_func=FEC_Util.name_func)

if __name__ == "__main__":
    run()
