# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,pickle,os,string

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities
from GeneralUtil.python import CheckpointUtilities

from Research.Personal.EventDetection.Util import Offline,Plotting

def get_best_metrics(data_file):
    learners = CheckpointUtilities.lazy_load(data_file)
    metrics = [Offline.best_metric_from_learner(l) for l in learners]
    return metrics

def precison_recall_plot(metric):
    precision = metric.precision()
    recall = metric.recall()
    f1_score = 2*precision*recall/(precision+recall)
    data = [precision,recall,f1_score]
    bins = np.linspace(0,1)
    colors = ['m','b','r']
    labels = ["Precision","Recall","F$_1$-Score"]
    ind = np.arange(len(data)) 
    style_precision = dict(color='k',alpha=0.3,hatch="//",label="Precision")
    style_recall = dict(color='r',alpha=0.3,label="Recall")
    width = 0.7
    rects = plt.bar(ind,data,color=colors,width=width,alpha=0.6)
    ax = plt.gca()
    ax.set_xticks(ind)
    ax.set_xticklabels(labels)
    y_func=lambda i,r: "{:.3f}".format(r.get_height()/2)
    PlotUtilities.autolabel(rects,y_func=y_func)
    PlotUtilities.lazyLabel("Metric",
                            "Value","",
                            frameon=True)
    plt.xlim([-width,len(labels)-width/2])
    plt.ylim([0,1])
    

def run(in_base="../FigurePerformance_CS/"):
    """
    
    """
    out_base = "./"
    data_file = "../FigurePerformance_CS/data/Scores.pkl"
    force=False
    cache_file = out_base + "cache.pkl"
    metrics = CheckpointUtilities.getCheckpoint(cache_file,get_best_metrics,
                                                force,data_file)
    loc_left = (-0.10,1.1)
    loc_top = (-0.12,1.05)
    loc_lower = (-0.12,0.95)
    locs = [loc_top,loc_left,loc_left,loc_lower,loc_lower]
    titles = Plotting.algorithm_title_dict()
    colors = Plotting.algorithm_colors()
    names = [titles[m.name.lower()] for m in metrics]
    # make the precision and recall plots
    fig = PlotUtilities.figure((12,4))
    for i,m in enumerate(metrics):
        ax = plt.subplot(1,3,(i+1))
        precison_recall_plot(m)
        PlotUtilities.title(names[i])
        ax.tick_params(axis=u'x', which=u'both',length=0)
    PlotUtilities.savefig(fig,"./f_score.pdf")
    for i,m in enumerate(metrics):
        name = names[i]
        safe_name = name.replace(" ","")
        color_pred =  colors[i]
        print("The best parameter for {:s} was {:.4g}".\
              format(name,m.x_values[m.best_param_idx]))
        distance_histogram= \
            Plotting.event_error_kwargs(m,color_pred=color_pred,
                                        q=Offline._def_q())
        true,pred = m.true,m.pred
        # make the 'just the distance' figures
        fig = PlotUtilities.figure((10,5))
        Plotting.histogram_event_distribution(use_q_number=True,
                                              **distance_histogram)
        final_out_dist = "{:s}{:s}_dist.pdf".format(out_base,safe_name)
        plt.ylim([0.5,max(plt.ylim())*3])
        PlotUtilities.savefig(fig,final_out_dist)
        # make the rupture spectrum figure
        fig = PlotUtilities.figure((14,7.5))
        final_out_rupture = "{:s}{:s}_rupture.pdf".format(out_base,safe_name)
        Plotting.rupture_plot(true,pred,fig=fig,
                              color_pred=color_pred)
        PlotUtilities.savefig(fig,final_out_rupture)
        fig = PlotUtilities.figure((16,8))
        # plot the metric plot
        Plotting.rupture_plot(true,pred,use_legend=True,
                              distance_histogram=distance_histogram,
                              fig=fig,color_pred=color_pred)
        final_out_path = "{:s}{:s}.pdf".format(out_base,safe_name)
        PlotUtilities.label_tom(fig,loc=locs,fontsize=18)
        plt.suptitle(name,fontsize=25,y=0.95,color=colors[i],alpha=0.7)
        PlotUtilities.savefig(fig,final_out_path,
                              subplots_adjust=dict(wspace=0.2,hspace=0.1,
                                                   left=0.05,top=0.85))



if __name__ == "__main__":
    run()
