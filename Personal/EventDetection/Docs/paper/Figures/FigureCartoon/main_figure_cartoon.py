# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os

import svgutils.compose as sc
sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection.Util import Analysis 
# /!\ note the 'SVG' function also in svgutils.compose
from IPython.display import SVG 

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
    fig = PlotUtilities.figure((12,4))
    for i,c in enumerate(cases):
        plt.subplot(1,len(cases),(i+1))
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
        PlotUtilities.tick_axis_number(num_x_major=5,num_y_major=3)
    out_tmp = "FigureCartoon{:d}.svg".format(i)
    out_names.append(out_tmp)
    w_space = 0.5
    PlotUtilities.savefig(fig,out_tmp,
                          subplots_adjust=dict(wspace=w_space,hspace=0))
    """
    see: 
    stackoverflow.com/questions/31452451/importing-an-svg-file-a-matplotlib-figure
    """
    tip_base = base  + "cartoon/2017-2-event-detection/" + \
               "SurfaceChemistry Dig10p3_combined_no_extra.svg"
    cartoon_files = [tip_base]
    tip_panels = [sc.Panel(sc.SVG(file_path)).scale(2.2)
                  for i,file_path in enumerate(cartoon_files)]
    data_panels = [sc.Panel(sc.SVG(f)) for f in out_names]
    all_panels = tip_panels + data_panels
    sc.Figure("41cm", "20cm", 
              *(tip_panels + data_panels)
    ).tile(1,2).save("final.svg")
    # remove all the intermediate svg files
    for f in out_names:
        os.remove(f) 

if __name__ == "__main__":
    run()
