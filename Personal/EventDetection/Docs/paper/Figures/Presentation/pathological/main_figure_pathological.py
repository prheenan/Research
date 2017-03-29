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

def plot_arrows(split_fec):
    x_plot,y_plot = Plotting.plot_format(split_fec.retract)
    marker_size = 30
    fudge_pN = 25
    event_idx = split_fec.get_retract_event_starts()
    Plotting.plot_arrows_above_events(event_idx,x_plot,y_plot,fudge_pN,
                                      markersize=marker_size,hatch='////')

def plot_retract_with_events(split_fec,fudge_pN=15,marker_size=30,
                             n_filter_points=1000):
    x_plot,y_plot = Plotting.plot_format(split_fec.retract)
    x_filtered,y_filtered = \
        FEC_Plot._fec_base_plot(x_plot,y_plot,n_filter_points=n_filter_points)
    plot_arrows(split_fec)
    
def plot_events_by_colors(split_fec,predicted_event_idx,n_points_filter=1000,
                          plot_filtered=False):
    start_idx = split_fec.get_retract_event_starts()
    x_plot,y_plot = Plotting.plot_format(split_fec.retract)
    retract_filtered = FEC_Util.GetFilteredForce(split_fec.retract,
                                                 n_points_filter)
    x_filtered,y_filtered = Plotting.plot_format(retract_filtered)
    fudge_negative = -30
    before =  [0] + list(predicted_event_idx)
    after =  list(predicted_event_idx) + [None]
    slices = [slice(i,f,1) for i,f in zip(before,after)]
    color_before = 'r'
    color_after = 'b'
    for i,(before_tmp,after_tmp) in enumerate(zip(slices[:-1],slices[1:])):
        dict_common = dict(color_before=color_before,
                           color_after=color_after,
                           before_slice=before_tmp,
                           after_slice=after_tmp)
        Plotting.before_and_after(x_plot,y_plot,style=dict(alpha=0.3),
                                  **dict_common)
        if (plot_filtered):
            Plotting.before_and_after(x_filtered,y_filtered,
                                      style=dict(alpha=1.0),**dict_common)
        tmp = color_before
        color_before = color_after
        color_after = tmp
    plot_arrows(split_fec)

def label_and_save(fig,name,num=0):
    PlotUtilities.lazyLabel("Time (s)","Force","")
    PlotUtilities.savefig(fig,name + "{:d}.pdf".format(num))
    return num+1

def run(base="./"):
    """
    
    """
    name = "examples.pdf"
    data_base = base + "data/"
    file_names = ["fast_unfolding","low-snr"]
    kw = dict(cache_directory=data_base,force=False)
    file_paths = [data_base + f +".csv" for f in file_names]
    cases = [read_and_cache_file(f,**kw) for f in file_paths]
    for c in cases:
        split_fec = Analysis.zero_and_split_force_extension_curve(c)
        name = c.Meta.Name
        num = 0
        retract_idx = split_fec.get_predicted_retract_surface_index()
        split_fec.retract.Force -= \
            np.percentile(split_fec.retract.Force[retract_idx:],25)
        predicted_event_idx = Detector.predict(c,threshold=1e-2)
        fig = PlotUtilities.figure((8,4))
        plot_retract_with_events(split_fec)
        num = label_and_save(fig=fig,name=name,num=num)
        fig = PlotUtilities.figure((8,4))
        plot_events_by_colors(split_fec,predicted_event_idx)
        num = label_and_save(fig=fig,name=name,num=num)
        fig = PlotUtilities.figure((8,4))
        plot_events_by_colors(split_fec,predicted_event_idx,plot_filtered=True)
        num = label_and_save(fig=fig,name=name,num=num)


if __name__ == "__main__":
    run()
