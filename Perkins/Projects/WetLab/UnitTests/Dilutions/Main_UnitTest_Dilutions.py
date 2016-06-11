# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../..")
from Util import DilutionUtil

def run():
    """
    tests the dilution functionality
    """
    stocks = np.array([90,76,62,40,110,100,70,30])
    volumes = 32
    desired = 30
    obj = DilutionUtil.GetDilutionObjects(stocks,volumes,desired)
    # expected volume we need to add
    expectedVtoAdd = [64.0,49.0666666667,34.1333333333,10.6666666667,
                      85.3333333333,74.6666666667,42.6666666667,0.0]
    actual = [o.AddVol for o in obj]
    np.testing.assert_allclose(actual,expectedVtoAdd)

if __name__ == "__main__":
    run()
