# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
baseDir = "../../"
sys.path.append(baseDir)

import Util.KmerUtil as KmerUtil
import Util.AlignUtil as AlignUtil
import Util.MeltingTemperatureUtil as melt

import PrimerDesign.OverHangingPrimers.OverhangGeneration as \
    OverhangGeneration

def run():
    """
    (1) gets a possible 'alphabet' of kmers (k=6) not appearing in the plasmid
    (2) converts the possible kmers to a 'psuedo plasmid', comprising sequences
        not in the plasmid
    (3) Runs emboss eprimer32, filters primers with melting temperatures we want
    (4) Filters by self-alignment score, using an algorithm to sort similar to 
    IDT
    
    In the end, the primers we find
    (1) Have no kmers in common with the plasmid
    (2) Have melting temperatures we want
    (3) Have a GC clamp, GC content maximized
    (4) Have low self-dimerization scores
    """
    inputFile = baseDir + "PlasmidData/Plasmids/mp13_plasmid_plasmid_seq.txt"
    # by experimentation, we know that all kmers of up to length five,
    # so start with 6 (could find this very easily...)
    kmerLen=6
    optTm=30
    primerLen =12
    mInfoObj = OverhangGeneration.GetBestOverhangs(inputFile,primerLen,optTm)
    mPlasmid = mInfoObj.Plasmid
    # Read in the plasmid file
    print("Best {:d}-mer is {:s}".format(3,KmerUtil.BestKmer(mPlasmid,3)))
    for seq in mInfoObj.BestPrimers:
        print("{:s} has dimer score {:.1f}, MT {:.2f}".\
              format(seq,seq.SelfAlignment,seq.temp))
    # use the actual primers
    forward = "GTGGTCCTAGTG"
    mPrimers = [forward,KmerUtil.ReverseComplement(forward)]
    print("Known primers...")
    for p in mPrimers:
        alignments = AlignUtil.AlignSelfWithReverseComplement(p)
        score = max(s.score for s in alignments)
        print("{:s} has dimer score {:.1f}".format(p,score))
        
if __name__ == "__main__":
    run()
