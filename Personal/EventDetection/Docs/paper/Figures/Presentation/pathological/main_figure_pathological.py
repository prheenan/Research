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

def plot_format(time_sep_force):
    x_plot = time_sep_force.Separation * 1e9
    y_plot = time_sep_force.Force * 1e12
    return x_plot,y_plot

def plot_arrows(split_fec):
    x_plot,y_plot = plot_format(split_fec.retract)
    marker_size = 30
    fudge_pN = 30
    event_idx = split_fec.get_retract_event_starts()
    Plotting.plot_arrows_above_events(event_idx,x_plot,y_plot,fudge_pN,
                                      markersize=marker_size,hatch='////')

def plot_retract_with_events(split_fec,fudge_pN=15,marker_size=30,
                             n_filter_points=1000):
    x_plot,y_plot = plot_format(split_fec.retract)
    x_filtered,y_filtered = \
        FEC_Plot._fec_base_plot(x_plot,y_plot,n_filter_points=n_filter_points)
    plot_arrows(split_fec)
    
def plot_events_by_colors(split_fec,predicted_event_idx):
    start_idx = split_fec.get_retract_event_starts()
    x_plot,y_plot = plot_format(split_fec.retract)
    fudge_negative = -30
    before =  [0] + list(predicted_event_idx)
    after =  list(predicted_event_idx) + [None]
    slices = [slice(i,f,1) for i,f in zip(before,after)]
    color_before = 'r'
    color_after = 'b'
    for i,(before_tmp,after_tmp) in enumerate(zip(slices[:-1],slices[1:])):
        Plotting.before_and_after(x_plot,y_plot,before_tmp,after_tmp,
                                  color_before=color_before,
                                  color_after=color_after,
                                  style=dict(alpha=0.5))
    tmp = color_before
    color_before = color_after
    color_after = tmp


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
        retract_idx = split_fec.get_predicted_retract_surface_index()
        split_fec.retract.Force -= \
            np.percentile(split_fec.retract.Force[retract_idx:],25)
        predicted_event_idx = Detector.predict(c,threshold=1e-2)
        fig = PlotUtilities.figure((8,4))
        plot_retract_with_events(split_fec)
        plot_events_by_colors(split_fec,predicted_event_idx)
        PlotUtilities.lazyLabel("Separation (nm)","Force","")
        PlotUtilities.savefig(fig,c.Meta.Name + ".pdf","_pres.pdf")

if __name__ == "__main__":
    run()
