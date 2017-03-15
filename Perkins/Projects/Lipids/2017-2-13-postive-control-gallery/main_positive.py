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
            "PositiveControls/SingleAttachments/"
    out = base + "RepresentativeGallery/"
    data_base = base + "Excessive_Concentrations/"
    threshold_meters = 10e-9
    adhesion_meters = 5e-9
    files,raw_data = FEC_Util.read_and_cache_pxp(data_base,force=False)
    proccess_opts = dict(zero_separation_at_zero_force=True)
    processed = [FEC_Util.SplitAndProcess(r,**proccess_opts) for r in raw_data]
    # make individual plots 
    retr_filtered = []
    proc_filtered = []
    for i,(appr,retr) in enumerate(processed):
        adhesion_idx = np.where(retr.Separation < adhesion_meters)[0]
        after_adhesion_idx = adhesion_idx[-1] + 1
        # only look where the retract max makes sense XXX should filter pxp
        arg_max = np.argmax(retr.Force[after_adhesion_idx:])
        if (retr.Separation[after_adhesion_idx:][arg_max] > threshold_meters):
            continue
        # POST: separation is some kind of sensible (less than threshold_meters)
        fig = PlotUtilities.figure()
        FEC_Plot.FEC_AlreadySplit(appr,retr,NFilterPoints=250)
        # get some sensible limit in nm/pN for x/y
        min_x = plt.xlim()[0]
        max_y = plt.ylim()[-1]
        max_x = min_x + 40
        plt.xlim([-5,max_x])
        plt.ylim([-30,max_y])
        PlotUtilities.savefig(fig,out+ "FEC_{:d}.pdf".format(i))
        # filter the retract to just the bits we care about
        idx_max = np.where(retr.Separation*1e9 >= max_x)[0][0]
        retr_sliced = FEC_Util.MakeTimeSepForceFromSlice(retr,slice(0,idx_max,1))
        retr_filtered.append(retr_sliced)
        proc_filtered.append(processed[i])
    # make a 2-D histogram of the fitered ones 
    fig = PlotUtilities.figure()
    FEC_Plot.heat_map_fec(retr_filtered)
    plt.ylim(-30,plt.ylim()[1])
    PlotUtilities.savefig(fig,out + "histogram.pdf")
    # make a gallery of the fitered ones 
    inches_per_plot = 4.5
    n_gallery = 12
    gallery_curves = proc_filtered[:n_gallery]
    n_rows,n_cols = FEC_Plot._n_rows_and_cols(gallery_curves)
    fig_size = (n_cols*inches_per_plot,n_rows*inches_per_plot)
    fig = PlotUtilities.figure(figsize=(fig_size))
    ylim_pN = [-20,75]
    xlim_nm = [-10,50]
    FEC_Plot.gallery_fec(gallery_curves,xlim_nm,ylim_pN)
    plt.suptitle("Negative Control Gallery",y=1.2,fontsize=25)
    PlotUtilities.savefig(fig,base + "out.png",close=False)
    PlotUtilities.savefig(fig,"./out.png")    
if __name__ == "__main__":
    run()
