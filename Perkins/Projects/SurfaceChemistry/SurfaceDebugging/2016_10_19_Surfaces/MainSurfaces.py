# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import PlotUtilities,CheckpointUtilities

import matplotlib.gridspec as gridspec
class SurfaceImage:
    """
    Class for encapsulating (XXX just the height) of an image
    """
    def __init__(self,Example,height):
        """
        Args:

            Example: a single ProcessSingleWave.XX object
            height: the heights extracted from example (Example is pretty much 
            exclusively used for meta information)
        """
        self.height = height
        self.pixel_size_meters = Example.ImagePixelSize()
        self.NumRows = height.shape[0]
        self.range_meters = self.pixel_size_meters * self.NumRows
    def height_nm(self):
        """
        Returns the height as a 2-D array in nm
        """
        return self.height * 1e9
    def height_nm_rel(self,pct_considered_surface = 5):
        """
        returns the height, relative to the 'surface' (see args) in nm

        Args:
             pct_considered_surface: the lowest pct heights are consiered to be
             the absolute surface. between 0 and 100
        Returns:
             height_nm_rel, offset to the pct
        """
        height_nm = self.height_nm()
        MinV = np.percentile(height_nm,pct_considered_surface)
        height_nm_rel = height_nm - MinV
        return height_nm_rel

def ReadImage(InFile):
    """
    Reads a *single* image from the given pxp file

    Args:
         InFile: path
    Returns:
         tuple of <Full Wave Object, height array>
    """
    Waves = PxpLoader.LoadAllWavesFromPxp(InFile,
                                          ValidFunc=PxpLoader.IsValidImage)
    Example = Waves[0]
    Data = Example.DataY
    height = Data[:,:,0]
    return Example,height

def ReadImageAsObject(file_name):
    """
    Convenience wrapper; reads the file as an image
    
    Args:
        file_name: See ReadImage
    Returns:
        SurfaceImage object
    """
    Example,height = ReadImage(file_name)
    return SurfaceImage(Example,height)

def PlotImage(Image,label,**kwargs):
    """
    Plots SurfaceImage as a greyscale map
    
    Args:
         Image:  output of ReadImageAsObject
         label: what to label (title) this with
         **kwargs: passed to imshow
    """
    height_nm_relative = Image.height_nm_rel()
    range_microns = Image.range_meters * 1e6
    plt.imshow( height_nm_relative,extent=[0,range_microns,0,range_microns],
                cmap=plt.cm.Greys,aspect='auto',**kwargs)
    # remove the ticks
    PlotUtilities.lazyLabel(r"Microns",r"Microns",label)
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right='off', left='off')
    PlotUtilities.colorbar("Height (nm)")

def PlotImageDistribution(Image,pct=95,bins=300,PlotLines=True,AddSigmas=True,
                          label="Height\nDistribution",**kwargs):
    """
    Plots the distribution of image heights, relative to the surface

    Args:
        Image: output of ReadImageAsObject
        label: see PlotImage
        pct: where we want to draw the line on the distribution, after the
        cummulative probability is over this number [0,100]
        
        bins: number of bins for the histogram
        kwargs: passed to hist
    Returns:
        tuple of <n,bins,patches>, se matplotlib.pyplot.hist
    """
    height_nm_relative = Image.height_nm_rel()
    n,bins,patches = plt.hist(height_nm_relative.ravel(),bins=bins,linewidth=0,
                              alpha=0.3,label=label,**kwargs)
    height_nm_rel_encompassing_pct = np.percentile(height_nm_relative,pct)
    if (PlotLines):
        plt.axvline(height_nm_rel_encompassing_pct,linewidth=3,color='r',
                    linestyle="--",
                    label=("{:d}%<={:.1f}nm".\
                           format(int(pct),height_nm_rel_encompassing_pct)))
    plt.gca().set_yscale('log')
    plt.ylim([0.5,max(n)*2])
    PlotUtilities.lazyLabel("Height above surface(nm)",
                            "Count","Distribution of Heights")
    return n,bins,patches

def MakePlot(SurfaceImage,label,OutPath,**kwargs):
    """
    Makes a simple plot of the desired distribution
    """
    # make plots
    fig = PlotUtilities.figure(figsize=(10/1.5,16/1.5))
    ax = plt.subplot(2,1,1)
    PlotImage(SurfaceImage,label,**kwargs)
    ax = plt.subplot(2,1,2)
    PlotImageDistribution(SurfaceImage)
    PlotUtilities.savefig(fig,OutPath)

def MakeGridPlot(ImageInfo,Limits,Base,Name,figsize):
    # first, loop through and make the 'normal' plots, per distribution
    Images = []
    OutBase = Base + "out/"
    InBase = Base + "in/"
    for file_name,label,kwargs in ImageInfo:
        Image = CheckpointUtilities.getCheckpoint(OutBase +"c_"+label+".pkl",
                                                  ReadImageAsObject,False,
                                                  InBase + file_name)
        MakePlot(Image,label,OutPath=OutBase + label + ".pdf",**kwargs)
        Images.append(Image)
    # now we make a grid
    NumRows = 3
    NumCols = len(ImageInfo)
    gs = gridspec.GridSpec(NumRows, NumCols)
    bins = []
    fig = PlotUtilities.figure(figsize=figsize)
    colors = ['r','g','b','k']
    for i in range(NumCols):
        im  = Images[i]
        label,kwargs = ImageInfo[i][1:]
        # plot these next to each other
        ax = plt.subplot(gs[0,i])
        PlotImage(im,label,**kwargs)
        ax = plt.subplot(gs[1,i])
        bins_tmp = PlotImageDistribution(im,color=colors[i])
        plt.xlim(Limits)
        bins.append(bins_tmp)
    ax = plt.subplot(gs[2, :])
    for i,im in enumerate(Images):
        bins_tmp = PlotImageDistribution(im,label=ImageInfo[i][1],
                                         color=colors[i],
                                         PlotLines=False)
        plt.xlim(Limits)
    PlotUtilities.savefig(fig,OutBase+Name+".pdf")
            

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
    Base = "/Volumes/group/4Patrick/Reports/" + \
           "2016-10-20-tip-and-surface-height-distributions/"
    files_tips =[
        ["1_EtchedLong.pxp","Etched Long",tip_kwargs],
        ["9-221-2016-functiuonalized-9-19-tStrept-in-pbs-ph7.4-with-mini.pxp","'Good' T-Strept Mini",tip_kwargs],
        ["2016-10-14-mini-pbs-t-strept-tip-batch-10-10-2016_1%_attachment.pxp","'Bad'(?) T-Strept Mini",tip_kwargs] ]
    files_surfaces = [
        ["2016-10-6-mini-pbs-koh-clened-glass-batch-10-5-2016.pxp","KOH",surface_kwargs],
        ["2016-10-18-mini-pbs-PEG-azide-surface-batch-0p10x_peg_0p015mg_mL_10-17-2016_0%_attachment_with_proteins_possibly-tips-fault.pxp",r"PEG-Azide-$\frac{1}{10}$x",surface_kwargs],
        ["2016-10-14-mini-pbs-PEG-azide-surface-batch-10-10-2016_1%_attachment.pxp",
         "PEG-Azide-1x",surface_kwargs],
        ["2016-10-18-mini-pbs-PEG-azide-surface-batch-3x_peg_0p45mg_mL_10-17-2016_0%_attachment_with_proteins_possibly-tips-fault.pxp","PEG-Azide-3x",surface_kwargs],
        ]
    MakeGridPlot(files_tips,[-50,200],Base,"TipGrid.pdf",figsize=(16,16))
    MakeGridPlot(files_surfaces,[-5,30],Base,"SurfaceGrid.pdf",
                 figsize=(24,16))

if __name__ == "__main__":
    run()
