# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os,string

sys.path.append("../../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
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
    file_names = ["fast_unfolding","low-snr"]
    kw = dict(cache_directory=data_base,force=False)
    file_paths = [data_base + f +".csv" for f in file_names]
    cases = [read_and_cache_file(f,**kw) for f in file_paths]
    for c in cases:
        fig = PlotUtilities.figure((8,4))
        split_fec = Analysis.zero_and_split_force_extension_curve(c)
        retract_filtered = FEC_Util.GetFilteredForce(split_fec.retract,
                                                     1000)
        plot_x = split_fec.retract.Separation * 1e9
        plot_y = split_fec.retract.Force * 1e12
        fudge_pN = 15
        start_idx = split_fec.get_retract_event_starts()
        event_idx = start_idx
        fudge_negative = -30
        marker_size = 30
        Plotting.plot_arrows_above_events(event_idx,plot_x,plot_y,fudge_pN,
                                          markersize=marker_size,hatch='////')
        predicted_event_idx = Detector.predict(c,threshold=1e-2)
        before =  [0] + list(predicted_event_idx)
        after =  list(predicted_event_idx) + [None]
        slices = [slice(i,f,1) for i,f in zip(before,after)]
        color_before = 'r'
        color_after = 'b'
        for i,(before_tmp,after_tmp) in enumerate(zip(slices[:-1],slices[1:])):
            Plotting.before_and_after(plot_x,plot_y,before_tmp,after_tmp,
                                      color_before=color_before,
                                      color_after=color_after,
                                      style=dict(alpha=0.5))
            tmp = color_before
            color_before = color_after
            color_after = tmp
        PlotUtilities.lazyLabel("Time(ms)","Force","")
        PlotUtilities.savefig(fig,c.Meta.Name + ".pdf","_pres.pdf")

if __name__ == "__main__":
    run()
