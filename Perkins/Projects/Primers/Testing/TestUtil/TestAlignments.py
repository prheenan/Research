# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../")

import Util.AlignUtil as AlignUtil

def TestEbiScores():
    """
    Test local alignment scores to mp13 plasmid using EBI DNAFUll matrix.

    See: ebi.ac.uk/Tools/psa/emboss_water/nucleotide.html
    """
    
    SeqAndScores = [
        # Following obtained 2016-8-18, using 'Matrix: EDNAFULL'
        ['GTACGACTAGGACCACTT',38.0],
        ['TACGACTAGGACCACTCT',40.5],
        ['GTACGACTAGGACCACTC',38],
        ['GTAGGCCTAGGACCACTT',38.5]
    ]
    
def TestReverseComplementAlignments():
    """
    Using data fed into http://www.idtdna.com/calc/analyzer,
    and using the self-dimer tool, gets the self-dimer scores for various
    sequences
    """
    seqAndScore = [
        # 3/<??>/2016
        ['TAGGACCACTCG',2],
        ['ACGACTAGGACC',4],
        ['ACCACTCGAGTG',10],
        ['GTGGTCCTAGTG',4],
        ['CACTAGGACCAC',4],
        ['ACCACTCGAGTG',10],
        # 4/26/2016
        ["CGGACCACTCTG",2],
        ["GGCAGAGTGGTCCTA",2],
        ["CGGGACCACTCT",2]
        #,
        #['GTGGTCCTAGGC',6] This actually breaks alignment -- problem with
        # biopython?  See stack overflow demo in Demos, also
#stackoverflow.com/questions/35874826/biopython-local-alignment-between-dna-sequences-doesnt-find-optimal-alignment
    ]
    for seq,score in seqAndScore:
        thisScore = AlignUtil.GetBestSelfDimerAlignmentScore(seq)
        assert np.allclose(thisScore,score)

def run():
    """
    Tests the alignment utilities
    """
    TestReverseComplementAlignments()

if __name__ == "__main__":
    run()
