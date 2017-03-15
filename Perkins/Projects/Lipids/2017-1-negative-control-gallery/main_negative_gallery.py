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
	

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    network = FEC_Util.default_data_root()
    base = network + "4Patrick/CuratedData/Lipids/DOPC/"+\
            "NegativeControls/Representative_Gallery/"
    _,raw_data = FEC_Util.read_and_cache_pxp(base,force=False)
    processed = [FEC_Util.SplitAndProcess(r) for r in raw_data]
    inches_per_plot = 4.5
    n_rows,n_cols = FEC_Plot._n_rows_and_cols(processed)
    fig_size = (n_cols*inches_per_plot,n_rows*inches_per_plot)
    fig = PlotUtilities.figure(figsize=(fig_size))
    ylim_pN = [-20,75]
    xlim_nm = [-10,50]
    FEC_Plot.gallery_fec(processed,xlim_nm,ylim_pN)
    plt.suptitle("Negative Control Gallery",y=1.2,fontsize=25)
    PlotUtilities.savefig(fig,base + "out.png",close=False)
    PlotUtilities.savefig(fig,"./out.png")

    # # plot each of the force extension curves in a separate subplot



if __name__ == "__main__":
    run()
