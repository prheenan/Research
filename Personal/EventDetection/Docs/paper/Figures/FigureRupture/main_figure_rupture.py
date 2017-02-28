# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection.Util import Analysis 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

def fmt(ax):
    ax.set_ylim([-30,60])
    ax.set_xlim([-15,500])

def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    out_fig = "out.png"
    example = read_and_cache_file(data_base + "rupture.csv",has_events=True,
                                  force=False,cache_directory=data_base)
    n_filter = 250
    kw = dict(cache_directory=data_base,force=False)
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
    style = dict(color='b')
    fig = PlotUtilities.figure((8,8))
    ax = plt.gca()
    # plot the force etc
    plt.plot(sep,force,alpha=0.3,**style)
    plt.plot(sep_filtered,force_filtered,**style)
    PlotUtilities.lazyLabel("Separation [nm]","Force (pN)","")
    plt.ylim([-30,60])
    plt.xlim([-15,300])
    fmt(ax)
    # plot the rupture
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    # zoom-factor: 2.5, location: upper-left
    left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
    bbox_to_anchor = [left,bottom,width,height]
    ax_insert = zoomed_inset_axes(ax, 4, loc=1)
    # zoom-factor: 2.5, location: upper-left
    ax_insert.set_xlim([150,200]) # apply the x-limits
    ax_insert.set_ylim([10,30]) # apply the y-limits
    plt.plot(sep,force,alpha=0.3,**style)
    plt.plot(sep_filtered,force_filtered,**style)
    mark_inset(ax, ax_insert, loc1=2, loc2=3, fc="none", ec="0.5")
    ax_insert.xaxis.set_visible('False')
    ax_insert.yaxis.set_visible('False')
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    PlotUtilities.savefig(fig,out_fig)
    

if __name__ == "__main__":
    run()
