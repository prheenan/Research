# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

import EmbossUtil

DataSubDir = "PlasmidData"
CircularSubDir = "CircularDNA"

def Get1607FAnd3520R(base="./"):
    """
    Args:
        base: where the primer ("PlasmidData") directory is located
    Returns:
        the 1607F and 3520R in 5'->3' as strings
    """
    dirSpace = "/"
    mPath = dirSpace.join([base,DataSubDir,CircularSubDir])
    mPath = mPath.replace(dirSpace+dirSpace,dirSpace) + dirSpace
    inFile1607 = "1607F_Normal.txt"
    inFile3520R = "3520R_Normal.txt"
    p1607F = EmbossUtil.ReadSimpleSeqFile(mPath + inFile1607)
    p3520R = EmbossUtil.ReadSimpleSeqFile(mPath + inFile3520R)
    return p1607F,p3520R

if __name__ == "__main__":
    run()
