# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


def run():
    """
    Runs the biopython shtick
    """
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    seq1 = "GTGGTCCTAGGC"
    seq2 = "GCCTAGGACCAC"
    match =1
    mismatch = -2
    gapopen = -2
    gapext = 0
    # http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
    # 'localms' takes <seq1,seq2, match.mismatch,open,extend>
    for a in pairwise2.align.localms(seq1,seq2,match,mismatch,gapopen,gapext):
        print(format_alignment(*a))
"""
This Prints:

GTGGTCCTAGGC----
      |||||
----GCCTAGGACCAC
  Score=5

But shouldnt we be able to get

GTGGTCCTAGGC----
     ||||||
----GCCTAGGACCAC
  Score=6 (ie: with the extra C-C aligned?
"""





if __name__ == "__main__":
    run()
