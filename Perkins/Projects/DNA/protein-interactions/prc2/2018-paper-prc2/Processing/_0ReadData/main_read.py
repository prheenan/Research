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
from Util import Processing

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    pass
    base_dir = Processing.raw_data_dir()
    cache_dir = Processing.cache_raw_images()
    images = ImageUtil.cache_images_in_directory(base_dir,cache_dir)    

if __name__ == "__main__":
    run()
