# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../")
from Util import KmerUtil,IdtUtil

def run():
    """
    Just going to write down the hairpin we want
    See SI, Table 2, Hairpin '30R50/4T' (n=10) of 
    Woodside, M. T. et al. 
    Nanomechanical measurements of the sequence-dependent folding landscapes of
    single nucleic acid hairpins. 
    Proc Natl Acad Sci U S A 103, 6190-6195 (2006).
    """
    Seq = "gagtcaacgtactgatcacgctggatcctaTTTTtaggatccagcgtgatcagtacgttgactc"
    Len = 25
    Forward = Seq[:Len]
    Reverse = Seq[-Len:]
    ReverseComplement = KmerUtil.ReverseComplement(Forward)
    assert ReverseComplement.lower() == Reverse.lower() , "Hairpin doesnt match"
    # POST: hairpin works, push it to IDT
    Name = "W_30R50_4T"
    SeqByChar = [l for l in Seq]
    Seqs =  [ SeqByChar,
              # add the DBCO to the 5' end, the Biotin to the 3' end
              [IdtUtil.Dbco5Prime()] + SeqByChar + [IdtUtil.Biotin3Prime()]]
    Options = dict(Scale=IdtUtil.Scales._100NM,
                   Purification=IdtUtil.Purifications.PAGE)
    # print off everything
    Seqs = IdtUtil.SequencesAndNamesToOrder(Seqs,[Name,Name+"_Labelled"],
                                            **Options)
    IdtUtil.PrintSequencesToFile(Seqs,Name + ".txt")

def AddDBCOAndBio():
    """
    Adds the DBCO and biotin strings under the following assumptions...
    (1) Biotin must be on end we pull from.
    End we pull from must be where 
    """
    
if __name__ == "__main__":
    run()
