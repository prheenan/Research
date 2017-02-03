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
    # write down Michaels hairpin sequence
    seq = "gag tca acg tac tga tca cgc tgg atc cta TTT Tta gga tcc agc gtg " +\
          "atc agt acg ttg act c"
    # add in thymines as spacrs
    seq = "TT" + seq + "TT"
    seq = seq.replace(" ","")
    labelled_seq = IdtUtil.AddDBCOAndBiotin(seq)
    Options = dict(Scale=IdtUtil.Scales._100NM,
                   Purification=IdtUtil.Purifications.PAGE)
    # print off everything
    Name = "W_48R50_4T_Lab"
    Seqs = IdtUtil.SequencesAndNamesToOrder([labelled_seq],
                                            [Name],
                                            **Options)
    IdtUtil.PrintSequencesToFile(Seqs,Name + ".txt")


if __name__ == "__main__":
    run()
