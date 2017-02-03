# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot

def read_pxp_in_dir(directory,**kwargs):
    files = GenUtilities.getAllFiles(directory,ext=".pxp")
    to_ret = []
    for f in files:
        raw_data = FEC_Util.ReadInData(f,**kwargs)
        to_ret.extend(raw_data)
    return to_ret

def read_and_cache_pxp(directory,cache_name =None,**kwargs):
    if (cache_name is None):
        cache_name = "./cache.pkl"
    return CheckpointUtilities.getCheckpoint(cache_name,read_pxp_in_dir,True,
                                             directory,**kwargs)

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    network = "/Volumes/group/"
    base = network + "4Patrick/CuratedData/Lipids/DOPC/"+\
            "NegativeControls/Representative_Gallery/"
    raw_data = read_and_cache_pxp(base)
    processed = [FEC_Util.SplitAndProcess(r) for r in raw_data]
    n_cols = 3
    n_rows = int(np.ceil(len(processed)/n_cols))
    inches_per_plot = 4.5
    fig_size = (n_cols*inches_per_plot,n_rows*inches_per_plot)
    fig = PlotUtilities.figure(figsize=(fig_size))
    ylim_pN = [-20,30]
    xlim_nm = [0,100]
    for i,r in enumerate(processed):
        plt.subplot(n_rows,n_cols,(i+1))
        appr,retr = r
        is_labelled = (i == 0)
        XLabel = "Stage Position (nm)" if is_labelled else ""
        YLabel = "Force (pN)"      if is_labelled else ""
        ApproachLabel = "Approach" if is_labelled else ""
        RetractLabel =  "Retract"  if is_labelled else ""
        PlotLabelOpts = dict(ApproachLabel=ApproachLabel,
                             RetractLabel=RetractLabel,
                             x_func = lambda x: x.ZSnsr)
        LegendOpts = dict(loc='upper left',frameon=True)
        FEC_Plot.FEC_AlreadySplit(appr,retr,XLabel=XLabel,YLabel=YLabel,
                                  LegendOpts=LegendOpts,
                                  PlotLabelOpts=PlotLabelOpts,NFilterPoints=100)
        plt.xlim(xlim_nm)
        plt.ylim(ylim_pN)
        ax = plt.gca()
        if (not is_labelled):
            ax.tick_params(labelbottom='off',labelleft='off')   
    plt.suptitle("Negative control image examples",y=1,fontsize=25)
    PlotUtilities.savefig(fig,base + "out.png")

    # # plot each of the force extension curves in a separate subplot



if __name__ == "__main__":
    run()
