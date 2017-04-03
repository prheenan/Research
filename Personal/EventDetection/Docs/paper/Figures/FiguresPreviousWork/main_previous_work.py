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
        plt.imshow(im,extent=extent,aspect="auto",origin=origin,
                   interpolation='bilinear')
        if (reversed_flag):
            plt.ylim
        if (i != n_subplots-1):
            PlotUtilities.no_x_label()
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
    out_dir = base
    label_sep = "Separation (nm)"
    label_none = ""
    label_force = "Force (pN)"
    def_fig = [6,4]
    # write down the figures, the x and y labels and x,y limits for each subplot
    figures = [ 
        # get the classificaiton figure
        [def_fig,
         [["Andreopoulos_2011_edit",label_sep,label_force,[0,70],[-10,275],
           False]]],
        # get the wavelet figure
        [def_fig,
         [["bentez_2017_edit_1",label_none,label_force,[0,550],[0,3000],True],
          ["bentez_2017_edit_2",label_sep,"Wavelet\nEnergy (au)",[0,550],
           [0,1],False]]],
        # get the contour length figure
        [def_fig,
         [["kuhn_2005_edit",label_sep,label_force,[-5,80],[-50,150],False]]],
        # get the water quality figure
        [def_fig,
         [["Perelman_2012_edit","Time (s)","Contaminants (mg/mL)",
           [0,4],[1.4,2.25],False]]],
        # get the stock price figure
        [def_fig,
         [["struzik_2002_price-vs-time_edit",
           "Time (seconds)","Stock price",
           [0,500],[116.4,118.2],False]]],
        # get the atmospheric figure
        [def_fig,
         [["turner_1994_edit",
           "Time (seconds)","Temperature\n(Celcius)",
           [0,500],[-1,2],False]]],
    ]
    for figsize,subplots in figures:
        fig = PlotUtilities.figure((figsize))
        plot_subplot(subplots,src)
        name =out_dir + str(subplots[0][0]) + ".pdf"
        PlotUtilities.savefig(fig,name)
    

if __name__ == "__main__":
    run()
