# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../")
from Util import KmerUtil,IdtUtil
import copy

def run():
    """
    "tacks" onto the woodside DNA construct to pad it to 100nt
    """
    DesiredLength = 68
    LoopLength = 4
    seq = "gag tca acg tac tga tca cgc tgg atc cta TTT Tta gga tcc agc gtg " +\
          "atc agt acg ttg act c"
    
    seq = "tt" + seq + "tt" 
    seq = seq.replace(" ","")
    original_seq = copy.deepcopy(seq)
    ExtraNeeded = DesiredLength-len(seq)
    assert (ExtraNeeded >= 0) , "Needed Positive number"
    assert (ExtraNeeded % 2) == 0 , "Must have even numbers to extend hairpin"
    length_each_end = int(ExtraNeeded/2)
    add_to_front = seq[:length_each_end]
    add_to_end = seq[-length_each_end:]
    seq = add_to_front + seq + add_to_end
    reverse_complement = KmerUtil.ReverseComplement(seq)
    stem_length = int((DesiredLength-LoopLength)/2)
    # do we now have the correct length?
    assert len(seq) == DesiredLength
    # make sure that the reverse complement of the hairpin is equal to the 
    # hairpin itself IE, we should have a match (excepting LOOP to LOOP')
    # 5' - F - LOOP - R  - 3'    (ie: forward sequence)
    # 3' - R'- LOOP'- F' - 5'    (ie: reverse complement)
    check = lambda s1,s2: s1.lower() == s2.lower()
    # is the reverse complement of the forward equal to the forward
    assert check(reverse_complement[:stem_length],
                  seq[:stem_length])
    # is the reverse complement of the reverse equal to the reverse
    assert check(reverse_complement[-stem_length:],
                 seq[-stem_length:])
    # the loops should *not* match
    assert not check(reverse_complement[stem_length:stem_length+LoopLength],
                     seq[stem_length:stem_length+LoopLength])
    labelled_seq = AddDBCOAndBiotin(seq)
    Options = dict(Scale=IdtUtil.Scales._100NM,
                   Purification=IdtUtil.Purifications.PAGE)
    # print off everything
    Name = "W_48R50_4T"
    Seqs = IdtUtil.SequencesAndNamesToOrder([seq,labelled_seq],
                                            [Name,Name+"_Labelled"],
                                            **Options)
    IdtUtil.PrintSequencesToFile(Seqs,Name + ".txt")


def AddDBCOAndBiotin(Sequence):
    SeqByChar = [l for l in Sequence]
    return [IdtUtil.Dbco5Prime()] + SeqByChar + [IdtUtil.Biotin3Prime()]


if __name__ == "__main__":
    run()
