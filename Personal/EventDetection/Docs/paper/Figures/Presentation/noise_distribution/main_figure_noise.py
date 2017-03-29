# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os,string

sys.path.append("../../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection.Util import Analysis,Plotting
# /!\ note the 'SVG' function also in svgutils.compose
import matplotlib.gridspec as gridspec
from Research.Personal.EventDetection._2SplineEventDetector import Detector

def run(base="./"):
    """
    
    """
    name = "examples.pdf"
    data_base = base + "data/"
    file_name = "fast_unfolding"
    kw = dict(cache_directory=data_base,force=False)
    file_path = data_base + file_name +".csv"
    fec = read_and_cache_file(file_path,**kw)
    split_fec = Analysis.zero_and_split_force_extension_curve(fec)
    x_plot,y_plot = Plotting.plot_format(split_fec.approach)
    interp = split_fec.approach_spline_interpolator()
    n_filter_points = 1000
    x_raw = split_fec.approach.Time
    y_raw = split_fec.approach.Force
    interp_raw = interp(x_raw)
    diff_raw = y_raw - interp_raw
    stdev = Analysis.local_stdev(diff_raw,n=split_fec.tau_num_points)
    f_plot_y = lambda y: y*1e12
    fig = PlotUtilities.figure((8,8))
    gs = gridspec.GridSpec(3,3)
    plt.subplot(gs[0,:])
    plt.plot(x_plot,y_plot,alpha=0.3,color='k')
    plt.plot(x_plot,f_plot_y(interp_raw),color='k',linewidth=3)
    PlotUtilities.lazyLabel("Time (s)","Force (pN)","")
    PlotUtilities.x_label_on_top()
    plt.subplot(gs[1,:])
    plt.plot(x_plot,f_plot_y(diff_raw),color='k',alpha=0.3)
    plt.plot(x_plot,f_plot_y(stdev),color='k',linewidth=3)
    PlotUtilities.lazyLabel("","Force (pN)","")
    PlotUtilities.no_x_label()
    plt.subplot(gs[2,0])
    PlotUtilities.no_x_anything()
    plt.subplot(gs[2,1])
    PlotUtilities.no_x_anything()
    plt.subplot(gs[2,2])
    PlotUtilities.no_x_anything()
    PlotUtilities.savefig(fig,"./out.png")


if __name__ == "__main__":
    run()
