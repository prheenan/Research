# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
baseDir = "../../"
sys.path.append(baseDir)

import Util.KmerUtil as KmerUtil
import Util.EmbossUtil as EmbossUtil
import Util.IdtUtil as IdtUtil



def run():
    """
    Reads in the normal annealing primers, adds a 3mer to each which appears
    the least in the plasmid and its reverse complement. 
    """
    primerFwdIn = baseDir +"PlasmidData/CircularDNA/1607F_Anneal_Overhang.txt"
    plasmidIn = baseDir +"PlasmidData/Plasmids/mp13_plasmid_plasmid_seq.txt"
    outDir = "./Output/"
    # Read in the overhanging primer file
    Forward1607F = EmbossUtil.ReadSimpleSeqFile(primerFwdIn)
    # Get the Reverse complement of the plasmid
    Reverse3520R = KmerUtil.ReverseComplement(Forward1607F)
    assert KmerUtil.ReverseComplement(Reverse3520R) == Forward1607F
    # read in the plasmid
    plasmid = EmbossUtil.ReadSimpleSeqFile(plasmidIn)
    # get a k=3mer to add (so we can order a competing product)
    # use the least occuring 3mer in the primer or the plasmid
    kCompete=3
    kmers = KmerUtil.LeastOccuringKmers(plasmid,kCompete)
    #POST: found a kmer not in the original primers
    leastOccuringKmer,NumberTimesOccuring = kmers[0]
    # we will make a 'competing' primer using the least occuring kmer
    # note that we use the same kmer for each
    competing1607F = leastOccuringKmer + Forward1607F
    competing3520R = leastOccuringKmer + Reverse3520R
    conv =lambda x: IdtUtil.IdtSeqStr(x)
    print("{:d}-mer {:s} occurs least ({:d}) times".\
          format(kCompete,conv(leastOccuringKmer),NumberTimesOccuring))
    print("1607F original: {:s}, Competing: {:s}".format(conv(Forward1607F),
                                                         conv(competing1607F)))
    print("3520R original: {:s}, Competing: {:s}".format(conv(Reverse3520R),
                                                         conv(competing3520R)))
    # save out competing sequences
    scale = IdtUtil.Scales._25NM
    purification = IdtUtil.Purifications.NONE
    mOrders = [
               IdtUtil.IdtOrder(ProductName="1607F_Competing_ANN",
                                Sequence=competing1607F,
                                Scale=scale,Purification=purification),
               IdtUtil.IdtOrder(ProductName="3520R_Competing_ANN",
                                Sequence=competing3520R,
                                Scale=scale,Purification=purification)
               ]
    IdtUtil.PrintSequencesToFile(mOrders,
                                 outDir + "CompetingCircularPrimers.txt")
                                
    

    
if __name__ == "__main__":
    run()
