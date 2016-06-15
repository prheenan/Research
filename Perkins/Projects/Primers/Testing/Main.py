# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../")
import TestUtil.MainTestUtil as MainUtil

def run():
    """
    Runs all the unit tests.
    """
    MainUtil.run()
    pass

if __name__ == "__main__":
    run()
