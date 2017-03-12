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
from Research.Personal.EventDetection._2SplineEventDetector import Detector

from Research.Personal.EventDetection.Util import Analysis 

def fmt(ax):
    ax.set_xlim([-0.1,2.5])
    ax.set_ylim([-50,40])

def run(base="./"):
    """
    
    """
    data_base = base + "data/"
    out_fig = "cartoon.png"
    example = read_and_cache_file(data_base + "rupture.csv",has_events=True,
                                  force=False,cache_directory=data_base)
    n_filter = 1000
    kw = dict(cache_directory=data_base,force=False)
    events = example.Events
    fec_split = Analysis.zero_and_split_force_extension_curve(example)
    event_idx_retract = fec_split.get_retract_event_centers()
    event_idx = event_idx_retract[0]
    zoom_factor = 20
    points_around = int(np.ceil(event_idx/zoom_factor))
    retract = fec_split.retract
    retract.Force -= np.median(retract.Force)
    # get everything in terms of ploting variables
    x_plot = lambda x: x - min(x)
    y_plot = lambda y: y * 1e12
    time = retract.Time
    x = x_plot(time)
    force = y_plot(retract.Force)
    interpolator = fec_split.retract_spline_interpolator()
    n_points = fec_split.tau_num_points
    force_filtered = y_plot(interpolator(time))
    stdev,_,_ = Analysis.stdevs_epsilon_sigma(force,force_filtered,
                                                        n_points)
    # get the epsilon and sigmas
    n = fec_split.approach.Force.size
    offset_approach =  n-fec_split.get_predicted_approach_surface_index()
    slice_fit_approach = slice(offset_approach,-offset_approach,1)
    epsilon,sigma = fec_split.calculate_epsilon_and_sigma(n_points=n_points,\
                                    slice_fit_approach=slice_fit_approach)
    epsilon_plot,sigma_plot = y_plot(epsilon),y_plot(sigma)
    #
    prob,_= Detector._no_event_probability(time,interp=interpolator,y=stdev,
                                           n_points=n_points,epsilon=epsilon,
                                           sigma=sigma)
    fig = PlotUtilities.figure((8,16))
    n_plots = 3
    plt.subplot(n_plots,1,1)
    plt.plot(x,force,color='k',alpha=0.3)
    plt.plot(x,force_filtered)
    PlotUtilities.lazyLabel("","Force [pN]","")
    plt.subplot(n_plots,1,2)
    plt.plot(x,stdev)
    plt.axhline(epsilon_plot,label="$\epsilon$")
    plt.axhline(epsilon_plot+sigma_plot,linestyle='--',
                label="$\epsilon \pm \sigma$")
    plt.axhline(epsilon_plot-sigma_plot,linestyle='--')
    PlotUtilities.lazyLabel("","R$_t$","")
    plt.subplot(n_plots,1,3)
    plt.plot(x,prob)
    PlotUtilities.lazyLabel("Time","R$_t$","")
    PlotUtilities.savefig(fig,out_fig)
    

if __name__ == "__main__":
    run()
