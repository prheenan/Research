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
    file_names = ["fast_unfolding"]
    kw = dict(cache_directory=data_base,force=False)
    file_paths = [data_base + f +".csv" for f in file_names]
    cases = [read_and_cache_file(f,**kw) for f in file_paths]
    fig = PlotUtilities.figure((8,4))
    for c in cases:
        split_fec = Plotting.plot_fec(c,use_events=False)
        plot_x = split_fec.retract.Separation * 1e9
        plot_y = split_fec.retract.Force * 1e12
        fudge_pN = 15
        start_idx = split_fec.get_retract_event_starts()
        event_idx = start_idx
        fudge_negative = -30
        marker_size = 30
        Plotting.plot_arrows_above_events(event_idx,plot_x,plot_y,fudge_pN,
                                          markersize=marker_size,hatch='////')
        predicted_event_idx = Detector.predict(c,threshold=1e-1)
        for i in predicted_event_idx:
            plt.axvline(plot_x[i])
        PlotUtilities.lazyLabel("Time(s)","Force","")
    PlotUtilities.savefig(fig,name.replace(".pdf","_pres.pdf"))

if __name__ == "__main__":
    run()
