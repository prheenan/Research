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
from Research.Personal.EventDetection.Util import Analysis 
# /!\ note the 'SVG' function also in svgutils.compose
import matplotlib.gridspec as gridspec

def plot_fec(example,color='r',n_filter=1000):
    fec_split = Analysis.zero_and_split_force_extension_curve(example)
    retract = fec_split.retract
    retract.Force -= np.median(retract.Force)
    retract_filtered = FEC_Util.GetFilteredForce(retract,n_filter)
    # get everything in terms of ploting variables
    x_plot = lambda x: x * 1e9
    y_plot = lambda y: y * 1e12
    sep = x_plot(retract.Separation)
    force = y_plot(retract.Force)
    sep_filtered = x_plot(retract_filtered.Separation)
    force_filtered = y_plot(retract_filtered.Force)
    style = dict(color=color)
    plt.plot(sep,force,alpha=0.3,**style)
    plt.plot(sep_filtered,force_filtered,**style)
    plt.tight_layout()

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
    styles = [dict(color='k'),
              dict(color='g'),
              dict(color='r')]
    fig_x_in = 16
    fig_y_in = 8
    im_path = base + "/cartoon/SurfaceChemistry Dig10p3_pmod-0{:d}.png"
    fig = PlotUtilities.figure((fig_x_in,fig_y_in))
    gs= gridspec.GridSpec(2,3)
    for i in range(3):
        plt.subplot(gs[0, i])
        image = plt.imread(im_path.format(i+1))
        plt.imshow(image,interpolation="nearest",aspect='auto',extent=None)
        ax = plt.gca()
        ax.axis('off')
    for i,c in enumerate(cases):
        plt.subplot(gs[1, i])
        style = styles[i]
        plot_fec(c,**style)
        not_first_plot = i != 0
        fmt(remove_y_labels=False,remove_x_labels=False)
        if (i == 0):
            y_label = "Force (pN)"
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
    letters =  string.lowercase[:n_categories]
    letters = [ ["({:s}{:d})".format(s,n+1) for s in letters]
                 for n in range(n_categories)]
    flat_letters = [v for list_of_v in letters for v in list_of_v]
    PlotUtilities.label_tom(fig,flat_letters,loc=(-1.19,1))
    PlotUtilities.savefig(fig,"cartoon.svg",
                          subplots_adjust=dict(wspace=0.4,hspace=0.1))

if __name__ == "__main__":
    run()
