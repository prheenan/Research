# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../")
from Util import KmerUtil,IdtUtil


def run():
    """
    """
    # generate a random stem with only 0.2% GC content
    desired_len = 100
    # add in GC_XX before stem
    str_v = "GC" * 8
    np.random.seed(42)
    seq_loop = "TTTT"
    len_spacers = 4
    len_loop = len(seq_loop)
    gc_start = "GC"
    n_start = len(gc_start)*2
    stem_size = int((desired_len - n_start - \
                     len(str_v) * 2-len_loop-len_spacers)/2)
    seq_stem = np.random.choice(a=["a","t","g","c"],p=[0.4,0.4,0.1,0.1],
                                size=stem_size)
    seq_stem = gc_start + "".join(seq_stem)
    seq = seq_stem + seq_loop + KmerUtil.ReverseComplement(seq_stem)
    stem_length = len(seq_stem)
    loop_length = len(seq) - 2 * stem_length
    # check the hairpin matches before we begin
    KmerUtil.assert_hairpin_matches(seq,stem_length)
    str_comp = KmerUtil.ReverseComplement(str_v)
    seq_loop = seq[stem_length:stem_length+loop_length]
    seq_gc_rich = seq[:stem_length] + str_v + seq_loop + str_comp + \
                  seq[-stem_length:]
    KmerUtil.assert_hairpin_matches(seq_gc_rich,stem_length+len(str_v))
    # make sure the sequences is the desired length
    # add in thymines as spacrs
    seq_full = "TT" + seq_gc_rich + "TT"
    # for a negative control, get the stem on the *3-prime* end (that is
    # where the biotin will live)
    n_stem_full = int((len(seq_full)-loop_length)/2)
    stem_one_side = KmerUtil.ReverseComplement(seq_full[-n_stem_full:])
    assert len(seq_full) == desired_len
    # print off everything
    Name = "W_100_GC18_4T"
    # get the unlablled versions
    Options = dict(Scale=IdtUtil.Scales._250NM)
    Seqs = IdtUtil.SequencesAndNamesToOrder([seq_full,stem_one_side],
                                            [Name,Name+"_3'Stem"],**Options)
    IdtUtil.PrintSequencesToFile(Seqs,Name + "_unlabelled.txt")
    # print the labelled versions
    labelled_seq = IdtUtil.AddDBCOAndBiotin(seq_full)
    Options = dict(Scale=IdtUtil.Scales._250NM,
                   Purification=IdtUtil.Purifications.PAGE)
    Name_Labelled = Name + "Lab"
    Seqs = IdtUtil.SequencesAndNamesToOrder([labelled_seq],
                                            [Name_Labelled],
                                            **Options)
    IdtUtil.PrintSequencesToFile(Seqs,Name + "_labelled.txt")


if __name__ == "__main__":
    run()
