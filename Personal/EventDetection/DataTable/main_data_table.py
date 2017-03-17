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



def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    positives_directory=  InputOutput.get_positives_directory()
    cache_directory = "../_1ReadDataToCache/cache/"
    categories = InputOutput.get_categories(positives_directory)
    force_read=False
    limit = 2
    categories = InputOutput.read_categories(categories,force_read,
                                             cache_directory,limit)
    for c in categores:
        print([d.Force.Size for d in c.data])

    pass

if __name__ == "__main__":
    run()
