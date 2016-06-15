# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../")
import Util.MeltingTemperatureUtil as melt

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    idt = lambda x : melt.GetIdtMeltingTemperature(x)
    # idt's specified error for DNA/DNA hybridization, see
    # http://www.idtdna.com/calc/analyzer
    abstol = 1.4
    # tested these sequences on IDTs website, with date as given
    seqsExpected = [
        # 12-mers (3/16/2016)
        ["AGA GTG GTC CTA",36.8],
        ["TAG GAC CAC TCT",36.8],
        ["TAG GAC CAC TCG",39.9],
        ["TAC GGA CCA CTC",40.4],
        ["CGG CGT GCG TCG",54.7],
        ["AGG AGA AGT GTC",36.2],
        ["AGT CTC CTT GTC",36.2],
        # 15-mers (3/18/2016)
        ["GTG TAG ACT GAA CTC",41.6],
        ["AGA GTG GTC CTA GAC",45.4],
        #  (4/26/2016)
        ["GCT ACG GAC ACT",41.5],
        ["CTA CGG ACC ACT CG",48.8],
        ["CGG ACC ACT CTG",43.1],
        ["GGC AGA GTG GTC CTA",49.8],
        ["GAC AGA GTG GTC CTA",45.8],
        ["CGG GAC CAC TCT",44.9]
    ]
    for seq,expected  in seqsExpected:
        actual = idt(seq)
        np.testing.assert_allclose(actual,expected,rtol=0,atol=abstol) 

if __name__ == "__main__":
    run()
