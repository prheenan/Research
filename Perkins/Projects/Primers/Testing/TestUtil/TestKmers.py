# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
sys.path.append("../../")

from collections import Counter
import Util.KmerUtil as KmerUtil

def TestKmerUtil():
    """
    Tests that the kmer utility works well
    """
    #
    # check getting the possible combinations
    #
    assert KmerUtil.AllPossibleKmers(1) == set(["A","C","G","T"])
    assert KmerUtil.AllPossibleKmers(2) == set(["AA","AC","AG","AT",
                                                "CA","CC","CG","CT",
                                                "GA","GC","GG","GT",\
                                                "TA","TC","TG","TT"])
    #
    # check getting complments and reverse complements
    #
    mStr = "ACTGCA"
    mStrComplement = "TGACGT"
    n = len(mStr)
    for i in range(n):
        # split the DNA into a 'forward' and 'back' portion
        fwdIdx = slice(0,i,1)
        backIdx = slice(i,None,1)
        fwd= mStr[:i]
        back = mStr[i:]
        # make the it getsthe complement right
        assert KmerUtil.Complement(fwd) == mStrComplement[fwdIdx]
        assert KmerUtil.Complement(back) == mStrComplement[backIdx]
        # make sure it gets the reverse complmenet right
        assert KmerUtil.ReverseComplement(fwd) == mStrComplement[fwdIdx][::-1]

    #
    #  Test that getting the kmers works fine
    #
    alpha = "ABCDEFGHIJKLMNOP"
    nAlpha = len(alpha)
    for i in range(nAlpha):
        kmers = [alpha[j:j+i] for j in range(nAlpha-i+1)]
        actual = KmerUtil.Kmers(alpha,i)
        actualSet = KmerUtil.KmerSet(alpha,i)
        assert kmers == actual
        assert set(kmers) == set(actual)
    #
    # test that getting the circularized kmers works fine
    # assumign that the kmers works fine
    #
    for i in range(1,len(mStr)):
        fullStr = mStr[-i:] + mStr + mStr[:i]
        nFull = len(fullStr)
        kmersFwd = [fullStr[j:j+i] for j in range(nFull-i)]
        revStr = KmerUtil.ReverseComplement(fullStr)
        kmersRev = [revStr[j:j+i] for j in range(nFull-i)]
        fullSet = set(kmersFwd) | set(kmersRev)
        # get all the possible ones
        allPossible = KmerUtil.AllPossibleKmers(i)
        # get the remainined kmers
        remaining = allPossible - fullSet
        actual = KmerUtil.CircularizedKmers(mStr,i)
        # check that just the circularized kmers works
        assert actual == fullSet
        # check that the 'remaining' kmers works
        actualNotAppearing = KmerUtil.GetKmersNotAppearing(mStr,i)
        assert actualNotAppearing == remaining
        # make sure we can get the correct number of counts
        allKmers = KmerUtil.Kmers(mStr,i)
        assert KmerUtil.GetKmerCounts(mStr,i) == Counter(allKmers)

def TestKmerPrimersValid():
    """
    Tests that we can correctly know if there is overlap between our primers
    """
    # want no more than 3 overlap
    kmerLen = 3
    mPlasmid = "ATCC"
    validKmers = ["GHI", # different alphabet
                  "ACG", # out of order
                  "AGG", # nope, skips
                  "1234", # no numbers in DNA...
                  "AT", # okay to have a 2-overlap
                  "GTT", # reverse, skips
              ]
    invalidKmers = [mPlasmid,
                    "ATC", # normal
                    "GAT", # reverse complement
                    "CAT", # crosses boundary (circular, forward)
                    "TGG", # crosses boundary (circular, reverse)
                    ]
    KmerUtil.AssertPrimersValid(mPlasmid,kmerLen,validKmers)
    # next do a negative case
    try:
        KmerUtil.AssertPrimersValid(mPlasmid,kmerLen,invalidKmers)
    except AssertionError as e:
        # as expected!
        return
    assert False , "Should have failed check"
        
def run():
    """
    Tests the kmer utility stuff
    """
    TestKmerUtil()
    TestKmerPrimersValid()

if __name__ == "__main__":
    run()
