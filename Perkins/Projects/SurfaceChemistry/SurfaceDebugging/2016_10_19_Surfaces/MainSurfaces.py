# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import PlotUtilities

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    pass
    Infile = "2016-10-6-mini-pbs-koh-clened-glass-batch-10-5-2016.pxp"
    Waves = PxpLoader.LoadAllWavesFromPxp(Infile,
                                          ValidFunc=PxpLoader.IsValidImage)
    Example = Waves[0]
    Data = Example.DataY
    height = Data[:,:,0]
    pixel_size_meters = Example.ImagePixelSize()
    NumRows = height.shape[0]
    range_meters = pixel_size_meters * NumRows
    range_microns = range_meters * 1e6
    height_nm = height * 1e9
    # look at the *relative* height
    pct_considered_surface = 5
    MinV = np.percentile(height_nm,pct_considered_surface)
    height_nm_rel = height_nm - MinV
    # make plots
    fig = PlotUtilities.figure(figsize=(5,8))
    ax = plt.subplot(2,1,1)
    plt.imshow( height_nm_rel,extent=[0,range_microns,0,range_microns],
                cmap=plt.cm.Greys,aspect='auto')
    # remove the ticks
    PlotUtilities.lazyLabel(r"Microns",r"Microns","Surface Image")
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right='off', left='off') 

    PlotUtilities.colorbar("Height (nm)")
    ax = plt.subplot(2,1,2)
    n,bins,patches = plt.hist(height_nm_rel.ravel(),bins=100,linewidth=0,
                              alpha=0.3,label="Height\nDistribution")
    pct = 95 
    height_nm_rel_encompassing_pct = np.percentile(height_nm_rel,pct)
    plt.axvline(height_nm_rel_encompassing_pct,linewidth=3,color='r',
                linestyle="--",
                label=("{:d}%<={:.1f}nm".\
                       format(pct,height_nm_rel_encompassing_pct)))
    ax.set_yscale('log')
    plt.ylim([0.5,max(n)*2])
    PlotUtilities.lazyLabel("Height above surface(nm)",
                            "Count","Distribution of Heights")
    PlotUtilities.savefig(fig,"out.png")

if __name__ == "__main__":
    run()
