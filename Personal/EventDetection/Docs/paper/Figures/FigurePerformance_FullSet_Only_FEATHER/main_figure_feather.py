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

def run(base="./"):
    """
    
    """
    out_base = base
    data_file = base + "data/Scores.pkl"
    force=False
    learners = CheckpointUtilities.lazy_load(data_file)
    for i,l in enumerate(learners):
        fig = PlotUtilities.figure((16,8))
        m = Offline.best_metric_from_learner(l)
        distance_histogram= Offline.event_error_kwargs(m)
        name,true,pred = m.name,m.true,m.pred
        # plot the metric plot
        Plotting.rupture_plot(true,pred,use_legend=True,
                              distance_histogram=distance_histogram,
                              fig=fig)
        final_out_path = "{:s}_{:d}.pdf".format(out_base,i)
        PlotUtilities.label_tom(fig,loc=(-0.25,1.1),fontsize=18)
        plt.suptitle("FEATHER",fontsize=25,y=0.95,fontcolor='b')
        PlotUtilities.savefig(fig,final_out_path,
                              subplots_adjust=dict(wspace=0.4,left=0.10,
                                                   top=0.8))




if __name__ == "__main__":
    run()
