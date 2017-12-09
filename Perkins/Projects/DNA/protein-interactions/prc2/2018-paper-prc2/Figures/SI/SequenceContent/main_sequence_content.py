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

def get_window_averaged_gc_content(f,n):
    s = SeqIO.read(f,"genbank").seq
    gc_content = np.zeros(len(s))
    chars = [c for c in s]
    gc = [1 if c == "C" or c == "G" else 0 for c in s]
    cumsum = np.cumsum(np.insert(gc,0,0))
    # gc[i] = mean(gc[i-n/2:i+n/2]))
    # in other words, gc[i] is the gc content *centered* at i 
    gc = (cumsum[n:] - cumsum[:-n])/n
    halfway = int(n/2)
    gc_content[halfway:-halfway+1] = gc
    # just have a constant boundary 
    gc_content[:halfway] = gc_content[halfway]
    gc_content[-halfway+1:] = gc_content[-halfway]
    return gc_content
    
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
    titles = ["12x601 Widom","GC-rich center","GC-poor center"]
    colors = ['r','g','b']
    gc_contents = [get_window_averaged_gc_content(f,n) for f in files]
    fig = PlotUtilities.figure((7,3.5))    
    for i,gc_content in enumerate(gc_contents):
        ax = plt.subplot(1,3,(i+1))
        plt.plot(gc_content,linewidth=0.75,color=colors[i])
        y_title = "<GC %> (Averaging {:d} bp)".format(n)
        PlotUtilities.lazyLabel("Position along sequence (bp)",y_title,titles[i])
        plt.ylim(0,1)
        if (i != 0):
            PlotUtilities.xlabel("")
            PlotUtilities.ylabel("")
            PlotUtilities.no_x_label(ax)
            PlotUtilities.no_y_label(ax)
    PlotUtilities.savefig(fig,"Sequence comparison")
    
    

if __name__ == "__main__":
    run()
