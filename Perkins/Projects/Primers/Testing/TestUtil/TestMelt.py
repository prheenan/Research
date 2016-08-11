# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../")
import Util.MeltingTemperatureUtil as melt

def TestMeltingTemperatures(Seqs,TemperatureFunction,abstol=1.4,rtol=0,
                            **kwargs):
    """
    Seq: list of tuples, like <Sequence, expected temperature in celcius>

    TemperatureFunction: function which takes in a sequence, gives back
    a temperature
     
    abstol: asbolute  temperature in celcius. Defaults to idt's specified error
    for DNA/DNA hybridization, see http://www.idtdna.com/calc/analyzer

    rtol: relative tolerance, defaults to zero
    kwargs: passed directly to assert_allclose
    """
    for i,(seq,expected)  in enumerate(Seqs):
        actual = TemperatureFunction(seq)
        print(i,seq,expected,actual,abs(expected-actual))
        np.testing.assert_allclose(actual,expected,atol=abstol,
                                   rtol=rtol,**kwargs)
        
def TestPCRMeltingTemperatures():
    idt = lambda x: melt.GetIdtMeltingTemperatureForPCR(x)
    SeqsExpected = [
        # 2016-8-10
        ["TAC GAC TAG GCC TAG AT",55],
        ["CTC CTA GTC GTA CGA CTA",56.1],
        ["AAG TGG TCC TAG TCG TAC",57.1],
        ]
    TestMeltingTemperatures(SeqsExpected,idt)
        
def TestBufferMeltingTemperatures():
    idt = lambda x : melt.GetIdtMeltingTemperature(x)
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
        ["CGG GAC CAC TCT",44.9],
        ["CGG GAC CAC TCT",44.9],
    ]
    TestMeltingTemperatures(seqsExpected,idt)


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    TestBufferMeltingTemperatures()
    TestPCRMeltingTemperatures()
    
if __name__ == "__main__":
    run()
