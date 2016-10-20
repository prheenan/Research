# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import PlotUtilities

class SurfaceImage:

    def __init__(self,Example,height):
        self.height = height
        self.pixel_size_meters = Example.ImagePixelSize()
        self.NumRows = height.shape[0]
        self.range_meters = self.pixel_size_meters * self.NumRows

def ReadImage(InFile):
    Waves = PxpLoader.LoadAllWavesFromPxp(InFile,
                                          ValidFunc=PxpLoader.IsValidImage)
    Example = Waves[0]
    Data = Example.DataY
    height = Data[:,:,0]
    return Example,height

def MakePlot(SurfaceImage,label,**kwargs):
    height_nm = SurfaceImage.height * 1e9
    # look at the *relative* height
    pct_considered_surface = 5
    MinV = np.percentile(height_nm,pct_considered_surface)
    height_nm_rel = height_nm - MinV
    # make plots
    fig = PlotUtilities.figure(figsize=(10/1.5,16/1.5))
    ax = plt.subplot(2,1,1)
    range_microns = SurfaceImage.range_meters * 1e6
    plt.imshow( height_nm_rel,extent=[0,range_microns,0,range_microns],
                cmap=plt.cm.Greys,aspect='auto',**kwargs)
    # remove the ticks
    PlotUtilities.lazyLabel(r"Microns",r"Microns",
                            "{:s} Surface Image".format(label))
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right='off', left='off') 
    PlotUtilities.colorbar("Height (nm)")
    ax = plt.subplot(2,1,2)
    n,bins,patches = plt.hist(height_nm_rel.ravel(),bins=1000,linewidth=0,
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
    PlotUtilities.savefig(fig,"{:s}.png".format(label))

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    surface_kwargs = dict(vmin=-1, vmax=4)
    tip_kwargs = dict(vmin=-10,vmax=18)
    in_files = [
        ["2016-10-6-mini-pbs-koh-clened-glass-batch-10-5-2016.pxp","KOH",surface_kwargs],
        ["2016-10-18-mini-pbs-PEG-azide-surface-batch-0p10x_peg_0p015mg_mL_10-17-2016_0%_attachment_with_proteins_possibly-tips-fault.pxp",r"PEG-Azide-$\frac{1}{10}$x",surface_kwargs],
        ["2016-10-14-mini-pbs-PEG-azide-surface-batch-10-10-2016_1%_attachment.pxp",
         "PEG-Azide-1x",surface_kwargs],
        ["2016-10-18-mini-pbs-PEG-azide-surface-batch-3x_peg_0p45mg_mL_10-17-2016_0%_attachment_with_proteins_possibly-tips-fault.pxp","PEG-Azide-3x",surface_kwargs],
        ["9-221-2016-functiuonalized-9-19-tStrept-in-pbs-ph7.4-with-mini.pxp","'Good' T-Strept Mini",
         tip_kwargs],
        ["2016-10-14-mini-pbs-t-strept-tip-batch-10-10-2016_1%_attachment.pxp","'Bad'(?) T-Strept Mini",
         tip_kwargs]]
    for file_name,label,kwargs in in_files:
        Example,height = ReadImage(file_name)
        Image = SurfaceImage(Example,height)
        MakePlot(Image,label,**kwargs)

if __name__ == "__main__":
    run()
