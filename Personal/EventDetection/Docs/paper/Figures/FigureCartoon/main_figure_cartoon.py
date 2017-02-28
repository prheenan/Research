# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

import svgutils.compose as sc
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection.Util import Analysis 
# /!\ note the 'SVG' function also in svgutils.compose
from IPython.display import SVG 

def plot_fec(example,color='r',n_filter=250):
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

def fmt(remove_x_labels=True,remove_y_labels=True):
    ax = plt.gca()
    if (remove_x_labels):
        ax.xaxis.set_ticklabels([])
    if (remove_y_labels):
        ax.yaxis.set_ticklabels([])
    PlotUtilities.lazyLabel("","","")
    plt.ylim([-30,60])
    plt.xlim([-15,500])

def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    kw = dict(cache_directory=data_base,force=False)
    multiple = read_and_cache_file(data_base + "multiple.csv",**kw)
    single = read_and_cache_file(data_base + "single.csv",**kw)
    no = read_and_cache_file(data_base + "no.csv",has_events=False,**kw)
    fig = PlotUtilities.figure((16,8))
    plt.subplot(1,3,1)
    plot_fec(no,color='k')
    fmt(remove_y_labels=False,remove_x_labels=False)
    PlotUtilities.ylabel("Force (pN)")
    PlotUtilities.xlabel("Separation (nm)")
    plt.subplot(1,3,2)
    plot_fec(single,color='g')
    fmt()
    plt.subplot(1,3,3)
    plot_fec(multiple,color='r')
    fmt()
    out_tmp = "FigureCartoon.svg"
    PlotUtilities.savefig(fig,out_tmp)
    # XXX need to fix this...
    # see: 
    #stackoverflow.com/questions/31452451/importing-an-svg-file-a-matplotlib-figure
    sc.Figure("16in", "16in", 
        sc.Panel(sc.SVG("./tip_attachments.svg").scale(3)),
        sc.Panel(sc.SVG(out_tmp)).move(0,600)
    ).save(out_tmp)
    

if __name__ == "__main__":
    run()
