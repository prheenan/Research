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
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    arr1 = np.loadtxt(f1,delimiter=",",skiprows=2)
    arr2 = np.loadtxt(f2,delimiter=",",skiprows=2)
    tol = 1e-6
    np.testing.assert_allclose(arr1,arr2,rtol=tol,atol=0,verbose=True)

if __name__ == "__main__":
    run()
