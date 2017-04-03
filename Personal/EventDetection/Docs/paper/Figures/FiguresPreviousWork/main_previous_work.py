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
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities


def plot_subplot(subplots,base):
    n_subplots = len(subplots)
    for i,(name,x,y,range_x,range_y,reversed_flag) in enumerate(subplots):
        plt.subplot(n_subplots,1,(i+1))
        in_file = base + name + ".png"
        im = plt.imread(in_file)
        extent = list(range_x) + list(range_y)
        if (not reversed_flag):
            origin= 'upper'
        else:
            origin = 'lower'
        plt.imshow(im,extent=extent,aspect="auto",origin=origin)
        if (reversed_flag):
            plt.ylim
        PlotUtilities.lazyLabel(x,y,"")

def run(base="./"):
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    src = base + "src/"
    label_sep = "Separation (nm)"
    label_none = ""
    label_force = "Force (pN)"
    # write down the figures, the x and y labels and x,y limits for each subplot
    figures = [ 
        # get the classificaiton figure
        [[8,8],
         [["Andreopoulos_2011_edit",label_sep,label_force,[0,70],[-10,275],
           False]]],
        # get the wavelet figure
        [[12,10],
         [["bentez_2017_edit_1",label_none,label_force,[0,550],[0,3000],True],
          ["bentez_2017_edit_2",label_sep,"Wavelet Energy (au)",[0,550],
           [0,1],False]]],
    ]
    for figsize,subplots in figures:
        fig = PlotUtilities.figure((figsize))
        plot_subplot(subplots,src)
        name = str(subplots[0][0]) + ".png"
        PlotUtilities.savefig(fig,name)
    

if __name__ == "__main__":
    run()
