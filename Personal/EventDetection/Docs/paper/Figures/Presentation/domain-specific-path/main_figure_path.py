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

def f_plot_y(y):
    return 1e12 * y

def make_plot(retract,pred_info,interp_raw,probability_idx,
              surface_idx,use_previous=True,use_surface_shading=True):
    # get the plotting information
    style_raw = dict(color='k',alpha=0.3)
    style_interp = dict(color='k',alpha=1)
    x_plot,y_plot = Plotting.plot_format(retract)
    min_x,max_x = min(x_plot),max(x_plot)
    fudge_x = abs(max_x-min_x) * 0.05
    xlim = np.array([min_x-fudge_x,max_x+fudge_x])
    probabilities = pred_info.probabilities
    probability_min = np.min([min(p) for p in probabilities])
    surface_plot = lambda label : \
        plt.axvspan(min(xlim),x_plot[surface_idx],alpha=0.2,color='r',
                    label=label)
    plt.subplot(2,1,1)
    plt.plot(x_plot,y_plot,**style_raw)
    plt.plot(x_plot,f_plot_y(interp_raw),**style_interp)
    if (use_surface_shading):
        surface_plot("Surface Contact")
    PlotUtilities.x_label_on_top()
    PlotUtilities.lazyLabel("Time (s)","Force (pN)","")
    plt.xlim(xlim)
    plt.subplot(2,1,2)
    if (use_previous and (probability_idx != 0)):
        plt.semilogy(x_plot,probabilities[probability_idx-1],color='k',
                     alpha=0.3,linestyle='--',label="Previous")
    plt.semilogy(x_plot,probabilities[probability_idx],label="Probability")
    plt.xlim(xlim)
    plt.ylim([probability_min/2,2])
    if (use_surface_shading):
        surface_plot(None)
    PlotUtilities.no_x_label()
    PlotUtilities.lazyLabel("","Probability","",loc="lower right")


def run(base="./"):
    """
    
    """
    name = "examples.pdf"
    data_base = base + "data/"
    file_name = "adhesions-2017-02-04-masters-data-650nm-dna-2.5ng_ul_1-15_dilution-25-hours-in-pbs-multiple-loading-rates-170203-1um_s_#1.pxpImage3516Concat"
    kw = dict(cache_directory=data_base,force=False)
    file_path = data_base + file_name +".csv"
    fec = read_and_cache_file(file_path,**kw)
    split_fec = Analysis.zero_and_split_force_extension_curve(fec)
    interp = split_fec.retract_spline_interpolator()
    retract = split_fec.retract
    interp_raw = interp(retract.Time)
    surface_idx = split_fec.get_predicted_retract_surface_index()
    # get the full prediction information
    _,pred_info = Detector._predict_full(fec)
    n_pred = 3
    base = "./pathway/"
    # make the initial figure, without shading
    fig = PlotUtilities.figure((6,10))
    make_plot(retract,pred_info,interp_raw,0,surface_idx,
              use_surface_shading=False)
    PlotUtilities.savefig(fig,base + "_initial.pdf")
    # make the intermediate figures
    for i in range(n_pred):
        fig = PlotUtilities.figure((6,10))
        make_plot(retract,pred_info,interp_raw,i,surface_idx)
        PlotUtilities.savefig(fig,base + "{:d}.pdf".format(i))
    # make the final figure
    fig = PlotUtilities.figure((6,10))
    make_plot(retract,pred_info,interp_raw,i,surface_idx,use_previous=False)
    PlotUtilities.savefig(fig,base + "{:d}.pdf".format(i+1))





if __name__ == "__main__":
    run()
