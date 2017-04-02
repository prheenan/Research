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

def plot_fec_cartoon(base,data_base,file_names,arrow_kwargs=dict()):
    kw = dict(cache_directory=data_base,force=False)
    file_paths = [data_base + f +".csv" for f in file_names]
    cases = [read_and_cache_file(f,**kw) for f in file_paths]
    n_cases = len(cases)
    out_names = []
    event_styles = Plotting._fec_event_colors
    styles = [dict(colors=event_styles,use_events=False),
              dict(colors=event_styles),
              dict(colors=event_styles)]
    fudge_pN = [10,12,25]
    im_path = base + "/cartoon/SurfaceChemistry Dig10p3_pmod-0{:d}.png"
    gs= gridspec.GridSpec(2,3)
    for i in range(3):
        plt.subplot(gs[0, i])
        image = plt.imread(im_path.format(i+1))
        plt.imshow(image,interpolation="bilinear",aspect='equal',extent=None)
        ax = plt.gca()
        ax.axis('off')
    for i,c in enumerate(cases):
        plt.subplot(gs[1, i])
        style = styles[i]
        fec_split = Plotting.plot_fec(c,**style)
        plt.xlim([-30,650])
        # decorate the plot to make it easier to read
        plot_x = fec_split.retract.Separation * 1e9
        plot_y = fec_split.retract.Force *1e12
        slices = fec_split.get_retract_event_slices()
        Plotting.top_bars(plot_x,plot_x,slices,colors=style['colors'])
        event_idx = [slice_v.stop for slice_v in slices]
        if (len(event_idx) > 0):
            # remove the last index (just te end of the FEC)
            event_idx = event_idx[:-1]
            fudge =fudge_pN[i]
            Plotting.plot_arrows_above_events(event_idx,plot_x,plot_y,fudge,
                                              **arrow_kwargs)
        not_first_plot = i != 0
        fmt(remove_y_labels=False,remove_x_labels=False)
        if (i == 0):
            y_label = r"Force (pN)"
            x_label = "Separation (nm)"
        else:
            y_label = ""
            x_label = ""
            ax = plt.gca()
            PlotUtilities.no_y_label(ax=ax)
        PlotUtilities.ylabel(y_label)
        PlotUtilities.xlabel(x_label)
        PlotUtilities.tick_axis_number(num_x_major=4)


def run(base="./"):
    """
    
    """
    name = "cartoon.pdf"
    data_base = base + "data/"
    file_names = ["no","single","multiple"]
    # save without the labels for the presentation
    subplots_adjust = dict(left=0.12,wspace=0.2,hspace=0.1)
    fig = PlotUtilities.figure((8,4))
    plot_fec_cartoon(base,data_base,file_names,
                     arrow_kwargs=dict(markersize=8))
    PlotUtilities.savefig(fig,name.replace(".pdf","_pres.pdf"),
                          subplots_adjust=subplots_adjust)
    # save with the labels for the presentation
    fig = PlotUtilities.figure((16,8))
    subplots_adjust = dict(left=0.08,wspace=0.2,hspace=0.1)
    plot_fec_cartoon(base,data_base,file_names)
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
    PlotUtilities.savefig(fig,name,subplots_adjust=subplots_adjust)

if __name__ == "__main__":
    run()
