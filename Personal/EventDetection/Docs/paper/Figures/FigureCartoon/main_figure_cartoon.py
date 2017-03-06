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
    for i,c in enumerate(cases):
        fig = PlotUtilities.figure((16,8))
        plt.subplot(1,3,1)
        style = styles[i]
        plot_fec(c,**style)
        not_first_plot = i != 0
        fmt(remove_y_labels=False,remove_x_labels=False)
        PlotUtilities.ylabel("Force (pN)")
        PlotUtilities.xlabel("Separation (nm)")
        plt.xlim([-30,650])
        out_tmp = "FigureCartoon{:d}.svg".format(i)
        out_names.append(out_tmp)
        PlotUtilities.savefig(fig,out_tmp)
    """
    see: 
    stackoverflow.com/questions/31452451/importing-an-svg-file-a-matplotlib-figure
    """
    tip_base = base  + "cartoon/2017-2-event-detection/" + \
               "SurfaceChemistry Dig10p3_pmod-0{:d}.svg"
    cartoon_files = [tip_base.format(i+1) for i in range(n_cases)]
    offset = 80
    delta = 20
    tip_panels = [sc.Panel(sc.SVG(file_path)).scale(2.5).move(offset*i+delta,0)
                  for i,file_path in enumerate(cartoon_files)]
    data_panels = [sc.Panel(sc.SVG(f)) for f in out_names]
    all_panels = tip_panels + data_panels
    sc.Figure("32cm", "32cm", 
              *(tip_panels + data_panels)
    ).tile(3, 2).save("final.svg")
    # remove all the files
    for f in out_names:
        os.remove(f) 
    

if __name__ == "__main__":
    run()
