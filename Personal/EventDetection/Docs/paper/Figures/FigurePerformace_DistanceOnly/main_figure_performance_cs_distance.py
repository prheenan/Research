# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle,os,string

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting,Analysis

def get_feather_run(data_file):
    return CheckpointUtilities.lazy_load(data_file)[0]

def make_distance_figure(l,data_file,fec_file):
    # get the multiple one to plot
    multiple_fec = CheckpointUtilities.lazy_load(fec_file)
    split_fec = Analysis.zero_and_split_force_extension_curve(multiple_fec)
    m = Offline.best_metric_from_learner(l)
    distance_histogram= Offline.event_error_kwargs(m)
    name,true,pred = m.name,m.true,m.pred
    event_starts = [i.start for i in split_fec.get_retract_event_slices()
                    if i.start is not 0]
    event_centers = split_fec.get_retract_event_centers()
    max_idx = max(event_centers)
    plot_x = split_fec.retract.Separation*1e9
    plot_y = split_fec.retract.Force*1e12
    plot_y -= np.median(plot_y)
    max_sep_nanometers = plot_x[max_idx]
    fudge = 1.5
    fudge_pN = 20
    plt.subplot(1,2,1)
    Plotting.plot_fec(multiple_fec)
    ylim = [-25,max(plt.ylim())]
    offset = min(plt.xlim())
    plt.xlim([-15,max_sep_nanometers*fudge])
    plt.ylim(ylim)
    pred_idx = int(min(event_starts)/2)
    # plot arrows from actual to predicted
    ax = plt.gca()
    head_length = 25
    fudge_arrow_pN = 5
    _arrowprops = dict(shrinkA=0,shrinkB=0.1,
                       width=4.0,alpha=0.3,edgecolor='k')
    arrowprops_true = dict(facecolor='g',hatch='//',**_arrowprops)
    arrowprops_pred = dict(facecolor='b',**_arrowprops)
    for start in event_starts:
        dy = 0
        dx = plot_x[pred_idx]-plot_x[start]
        xy = (plot_x[start],plot_y[start]+fudge_pN)
        xy_offset = (xy[0] + dx,xy[1] + dy)
        ax.annotate("", xy=xy_offset, xytext=xy,arrowprops=arrowprops_true)
    text_loc_x = np.mean([plot_x[min(event_starts)],plot_x[pred_idx]])
    text_loc_y = xy[1] + fudge_arrow_pN
    ax.text(text_loc_x,text_loc_y,s=r"d$_{\mathrm{t}\rightarrow\mathrm{p}}$",
            fontsize=20,color=arrowprops_true['facecolor'])
    # plot arrows from predicted to true
    y_pred = (max([(plot_y[pred_idx] + fudge_pN) for start in event_starts]) + \
              fudge_arrow_pN)
    closest_idx_true = min(event_starts)
    closest_x_true = plot_x[closest_idx_true]
    x_pred = plot_x[pred_idx]
    xy_pred = (x_pred,y_pred)
    xy_offset = (closest_x_true,y_pred)
    ax.annotate("", xy=xy_offset, xytext=xy_pred,arrowprops=arrowprops_pred)
    ax.text(x_pred*1.05,y_pred-fudge_arrow_pN*0.8,
            s=r"d$_{\mathrm{p}\rightarrow\mathrm{t}}$",
            fontsize=20,color=arrowprops_pred['facecolor'])
    plt.axvline(plot_x[pred_idx],linestyle='--',color='b',
                label="Predicted event")
    # plot the events as arrows
    Plotting.plot_arrows_above_events(event_starts,plot_x,plot_y,fudge_pN)
    PlotUtilities.lazyLabel("Separation (nm)","Force(pN)","",
                            frameon=True)
    plt.subplot(1,2,2)
    Plotting.histogram_event_distribution(**distance_histogram)
    ylim = plt.ylim()
    plt.ylim(min(ylim)/2,max(ylim)*2)


def run(base="./"):
    """
    
    """
    out_base = base
    data_base = base + "data/"
    data_file = data_base + "Scores.pkl"
    force=False
    cache_file = data_base + "cache.pkl"
    fec_file = data_base + "multiple.csv.pkl"
    name = "FEATHER"
    l = CheckpointUtilities.getCheckpoint(cache_file,get_feather_run,force,
                                          data_file)
    # make the rupture spectrum figure
    fig = PlotUtilities.figure((16,6))
    final_out_rupture = "{:s}{:s}_rupture.pdf".format(out_base,name)
    metric = Offline.best_metric_from_learner(l)
    x,name,true,pred = metric.x_values,metric.name,metric.true,metric.pred
    Plotting.rupture_plot(true,pred,fig=fig)
    PlotUtilities.savefig(fig,final_out_rupture)
    # make the distance histogram figure
    fig = PlotUtilities.figure((16,6))
    make_distance_figure(l,data_file,fec_file)
    final_out_hist = "{:s}{:s}_distances.pdf".format(out_base,
                                                     name.replace(" ","_"))
    PlotUtilities.label_tom(fig,loc=(-0.1,1.0),fontsize=18)
    PlotUtilities.savefig(fig,final_out_hist)



if __name__ == "__main__":
    run()
