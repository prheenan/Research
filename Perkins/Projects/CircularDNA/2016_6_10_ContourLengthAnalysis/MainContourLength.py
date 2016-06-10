# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../")


def run():
    """
    Runs contour length analysis
    """
    InFile = ".pxp"
    OutFile = ""
    Limit = 2
    # Read in the pxp (assume each 'name-group' with the same numerical
    # suffix represents a valid wave with a WLC of interest)

    # get just 'Limit' of the waves

    # get just the retract, normalize the force and separation to zeros.

    # fit an extensible worm-like chain to each, using grid search at a
    # 0.1nm scale for contour length and 0.001nm scale for persistence length
    # within 10% of expected (600-700nm,40-50 for persistence)
    # save this course-grained information

    # plot a histogram of the results

    
    
    

if __name__ == "__main__":
    run()
