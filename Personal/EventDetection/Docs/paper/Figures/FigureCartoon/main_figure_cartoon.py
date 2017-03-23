# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os,string

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection.Util import Analysis,Plotting
# /!\ note the 'SVG' function also in svgutils.compose
import matplotlib.gridspec as gridspec

def fmt(remove_x_labels=True,remove_y_labels=True):
    ax = plt.gca()
    if (remove_x_labels):
        ax.xaxis.set_ticklabels([])
    if (remove_y_labels):
        ax.yaxis.set_ticklabels([])
    PlotUtilities.lazyLabel("","","")
    plt.ylim([-30,60])
    plt.xlim([-30,700])

def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    kw = dict(cache_directory=data_base,force=False)
    file_names = ["no","single","multiple"]
    file_paths = [data_base + f +".csv" for f in file_names]
    cases = [read_and_cache_file(f,**kw) for f in file_paths]
    n_cases = len(cases)
    out_names = []
    styles = [dict(colors='r',use_events=False),
              dict(colors=['r','b']),
              dict(colors=['r','b','k'])]
    fig_x_in = 16
    fig_y_in = 8
    im_path = base + "/cartoon/SurfaceChemistry Dig10p3_pmod-0{:d}.png"
    fig = PlotUtilities.figure((fig_x_in,fig_y_in))
    gs= gridspec.GridSpec(2,3)
    for i in range(3):
        plt.subplot(gs[0, i])
        image = plt.imread(im_path.format(i+1))
        plt.imshow(image,interpolation="nearest",aspect='equal',extent=None)
        ax = plt.gca()
        ax.axis('off')
    for i,c in enumerate(cases):
        plt.subplot(gs[1, i])
        style = styles[i]
        Plotting.plot_fec(c,**style)
        not_first_plot = i != 0
        fmt(remove_y_labels=False,remove_x_labels=False)
        if (i == 0):
            y_label = r"Force (pN)"
            x_label = "Separation (nm)"
        else:
            y_label = ""
            x_label = ""
            ax = plt.gca()
            PlotUtilities.no_y_ticks(ax=ax)
        PlotUtilities.ylabel(y_label)
        PlotUtilities.xlabel(x_label)
        plt.xlim([-30,650])
        PlotUtilities.tick_axis_number(num_x_major=4)
    n_subplots = 2
    n_categories = len(file_names)
    letters =  string.uppercase[:n_categories]
    letters = ([r"{:s}".format(s) for s in letters] + \
               ["" for _ in range(n_categories)])
    bottom = (-0.25,1)
    top = (-0.60,1)
    loc = [top for i in range(n_categories)] +  \
          [bottom for i in range(n_categories)] 
    PlotUtilities.label_tom(fig,letters,loc=loc)
    PlotUtilities.savefig(fig,"cartoon.pdf",
                          subplots_adjust=dict(left=0.08,wspace=0.2,hspace=0.1))

if __name__ == "__main__":
    run()
