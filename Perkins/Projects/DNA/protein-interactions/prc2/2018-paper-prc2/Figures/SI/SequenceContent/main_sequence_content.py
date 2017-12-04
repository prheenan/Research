# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../../../../")
from GeneralUtil.python import PlotUtilities,GenUtilities
from Bio import SeqIO

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    data_base = "../../../Data/Sequences/"
    ext = "gb"
    files = GenUtilities.getAllFiles(data_base,ext)
    n = 50
    for f in files:
        s = SeqIO.read(f,"genbank").seq
        gc_content = np.zeros(len(s))
        chars = [c for c in s]
        gc = [1 if c == "C" or c == "G" else 0 for c in s]
        cumsum = np.cumsum(np.insert(gc,0,0))
        gc = (cumsum[n:] - cumsum[:-n])/n
        halfway = int(n/2)
        gc_content[halfway:-halfway+1] = gc
        gc_content[:halfway] = gc_content[halfway]
        gc_content[-halfway+1:] = gc_content[-halfway]
        average = np.mean(gc_content)
        fig = PlotUtilities.figure()
        plt.plot(gc_content)
        plt.axhline(average)
        PlotUtilities.lazyLabel("Position along sequence (bp)","GC fraction","")
        PlotUtilities.savefig(fig,f+".png")
    
    

if __name__ == "__main__":
    run()
