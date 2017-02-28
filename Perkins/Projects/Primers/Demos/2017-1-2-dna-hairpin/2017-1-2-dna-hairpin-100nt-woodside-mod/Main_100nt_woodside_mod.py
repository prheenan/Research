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
    stem_len = len(seq)-LoopLength
    KmerUtil.assert_hairpin_matches(seq,stem_length=stem_len)
    length_each_end = int(ExtraNeeded/2)
    add_to_front = seq[:length_each_end]
    add_to_end = seq[-length_each_end:]
    seq = add_to_front + seq + add_to_end
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
