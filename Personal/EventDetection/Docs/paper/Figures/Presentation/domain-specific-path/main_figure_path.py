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

def make_plot(retract,pred_info,interp_raw,probability_idx):
    # get the plotting information
    style_raw = dict(color='k',alpha=0.3)
    style_interp = dict(color='k',alpha=1)
    x_plot,y_plot = Plotting.plot_format(retract)
    xlim = np.array([min(x_plot),max(x_plot)])
    probabilities = pred_info.probabilities
    probability_min = np.min([min(p) for p in probabilities])
    plt.subplot(2,1,1)
    plt.plot(x_plot,y_plot,**style_raw)
    plt.plot(x_plot,f_plot_y(interp_raw),**style_interp)
    PlotUtilities.x_label_on_top()
    PlotUtilities.lazyLabel("Time (s)","Force (pN)","")
    plt.xlim(xlim)
    plt.subplot(2,1,2)
    for i in probability_idx:
        plt.semilogy(x_plot,probabilities[i])
    plt.xlim(xlim)
    plt.ylim([probability_min/2,2])
    PlotUtilities.no_x_label()
    PlotUtilities.lazyLabel("","Probability","")


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
    # get the full prediction information
    _,pred_info = Detector._predict_full(fec)
    n_pred = 3
    # get the limits
    for i in range(n_pred):
        fig = PlotUtilities.figure()
        make_plot(retract,pred_info,interp_raw,[i])
        PlotUtilities.savefig(fig,"./out{:d}.pdf".format(i))





if __name__ == "__main__":
    run()
