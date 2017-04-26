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
from Research.Personal.EventDetection.Util import Analysis
	
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
            "PositiveControls/SingleAttachments/"
    out = base + "RepresentativeGallery/"
    data_base = base + "10^4/"
    files,raw_data = FEC_Util.read_and_cache_pxp(data_base,force=False)
    processed = [Analysis.zero_and_split_force_extension_curve(r) 
                 for r in raw_data]
    # make individual plots 
    retr_filtered = []
    proc_filtered = []
    for i,split in enumerate(processed):
        # POST: separation is some kind of sensible (less than threshold_meters)
        fig = PlotUtilities.figure()
        appr,retr = split.approach,split.retract
        FEC_Plot.FEC_AlreadySplit(appr,retr,NFilterPoints=250)
        offset_meters = min(retr.Separation) + 15e-9
        appr.Separation -= offset_meters
        retr.Separation -= offset_meters
        # get some sensible limit in nm/pN for x/y
        min_x = plt.xlim()[0]
        max_y = plt.ylim()[-1]
        plt.ylim([-30,max_y])
        PlotUtilities.savefig(fig,out+ "FEC_{:d}.png".format(i))
        # filter the retract to just the bits we care about
        retr_filtered.append(retr)
        proc_filtered.append([appr,retr])
    # make a 2-D histogram of the fitered ones 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(retr_filtered,separation_max=60)
    plt.ylim(-30,plt.ylim()[1])
    PlotUtilities.savefig(fig,out + "histogram.png")
    # make a gallery of the fitered ones 
    inches_per_plot = 4.5
    n_gallery = 12
    gallery_curves = proc_filtered[:n_gallery]
    n_rows,n_cols = FEC_Plot._n_rows_and_cols(gallery_curves)
    fig_size = (n_cols*inches_per_plot,n_rows*inches_per_plot)
    fig = PlotUtilities.figure(figsize=(fig_size))
    ylim_pN = [-20,75]
    xlim_nm = [-30,75]
    FEC_Plot.gallery_fec(gallery_curves,xlim_nm,ylim_pN)
    plt.suptitle("Gallery",y=1.2,fontsize=25)
    PlotUtilities.savefig(fig,out + "gallery.png",close=False)
    PlotUtilities.savefig(fig,"./out.png")    
if __name__ == "__main__":
    run()
