# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
sys.path.append("../../")

import TestMelt
import TestAlignments
import TestKmers


def run():
    """
    Tests all the utility functions
    """
    TestKmers.run()
    TestAlignments.run()
    TestMelt.run()
    pass

if __name__ == "__main__":
    run()
