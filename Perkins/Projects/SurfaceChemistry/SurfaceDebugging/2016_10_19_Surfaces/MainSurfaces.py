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
         List of tuple of <Full Wave Object, height array>
    """
    Waves = PxpLoader.LoadAllWavesFromPxp(InFile,
                                          ValidFunc=PxpLoader.IsValidImage)
    # get all the images
    return [ (Example,Example.DataY[:,:,0]) for Example in Waves]

def ReadImageAsObject(file_name):
    """
    Convenience wrapper; reads the file as an image
    
    Args:
        file_name: See ReadImage
    Returns:
        List of SurfaceImage objects present in the file
    """
    Tuples = ReadImage(file_name)
    return [SurfaceImage(ex,height) for ex,height in Tuples]

def PlotImage(Image,**kwargs):
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
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right='off', left='off')

def PlotImageDistribution(Image,pct=95,bins=300,PlotLines=True,AddSigmas=True,
                          **kwargs):
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
                              edgecolor="none",alpha=0.3,
                              normed=False,**kwargs)
    height_nm_rel_encompassing_pct = np.percentile(height_nm_relative,pct)
    if (PlotLines):
        plt.axvline(height_nm_rel_encompassing_pct,linewidth=3,color='r',
                    linestyle="--",
                    label=("{:d}%<={:.1f}nm".\
                           format(int(pct),height_nm_rel_encompassing_pct)))
    plt.gca().set_yscale('log')
    plt.ylim([min(n)/2,max(n)*2])
    return n,bins,patches

def MakePlot(SurfaceImage,label,OutPath,**kwargs):
    """
    Makes a simple plot of the desired distribution
    """
    # make plots
    fig = PlotUtilities.figure(figsize=(10/1.5,16/1.5))
    ax = plt.subplot(2,1,1)
    PlotImage(SurfaceImage,**kwargs)
    PlotUtilities.lazyLabel(r"Microns",r"Microns",label)
    PlotUtilities.colorbar("Height (nm)")
    ax = plt.subplot(2,1,2)
    bins_tmp = PlotImageDistribution(SurfaceImage)
    PlotUtilities.lazyLabel("Height above surface(nm)",
                            "Count","Distribution of Heights")
    PlotUtilities.savefig(fig,OutPath)
    return bins_tmp

def MakeGridPlot(ImageInfo,Limits,Base,Name,figsize,Force=False,
                 ImageFunc=lambda m_list: m_list[0],hist_kwargs=dict()):
    """
    Given a list of images, files, and labels, makes a gridwise comparison

    Args:
        ImageInfo: Tuple of <FileName,Label,plotkwargs>
        Limits: x Limits for the histogram
        Base: Base directory (where in and out live
        figsize: how big to make the grid figure
        Force: if true, force re-reading
        ImageFunc: given a list of images from the file, selects one and 
        only one
    """
    # first, loop through and make the 'normal' plots, per distribution
    Images = []
    OutBase = Base + "out/"
    InBase = Base + "in/"
    max_n = 0
    for i,(file_name,label,kwargs) in enumerate(ImageInfo):
        List = CheckpointUtilities.getCheckpoint(OutBase +"c_"+label+".pkl",
                                                 ReadImageAsObject,Force,
                                                 InBase + file_name)
        Image = ImageFunc(List)
        n,bins,patches = \
            MakePlot(Image,label,OutPath=OutBase + str(i)+"_" +label + ".pdf",
                     **kwargs)
        max_n = max(max_n,max(n))
        Images.append(Image)
    # now we make a grid
    NumRows = 4
    NumCols = len(ImageInfo)
    gs = gridspec.GridSpec(NumRows, NumCols)
    bins = []
    fig = PlotUtilities.figure(figsize=figsize)
    colors = ['r','g','b','k','m','c','y']
    for i in range(NumCols):
        im  = Images[i]
        color_idx= i % len(colors)
        label,kwargs = ImageInfo[i][1:]
        # plot these next to each other
        ax = plt.subplot(gs[0,i])
        PlotImage(im,**kwargs)
        if (i ==0 ):
            PlotUtilities.lazyLabel(r"Microns",r"Microns",label)
            PlotUtilities.colorbar("Height (nm)")
        else:
            PlotUtilities.lazyLabel(r"","",label)
            PlotUtilities.colorbar("")
        ax = plt.subplot(gs[1,i])
        bins_tmp = PlotImageDistribution(im,color=colors[color_idx])
        plt.xlim(Limits)
        # set common y limits
        plt.ylim(0.5,2*max_n)
        bins.append(bins_tmp)
        # only add the x and y labels to the first plot, to de-clutter
        if (i == 0):
            PlotUtilities.lazyLabel("Height above surface(nm)",
                                    "Count","Distribution of Heights")
        else:
            # just get the axis formatting
            PlotUtilities.lazyLabel("","","")
    ax = plt.subplot(gs[2, :])
    labels = [i[1] for i in ImageInfo]
    for i,im in enumerate(Images):
        color_idx = i % len(colors)
        PlotImageDistribution(im,label=labels[i],color=colors[color_idx],
                              PlotLines=False,**hist_kwargs)
        plt.xlim(Limits)
        PlotUtilities.lazyLabel("Height above surface(nm)",
                                "Count","Distribution of Heights")
    ax = plt.subplot(gs[3,:])
    data = [i.height_nm_rel() for i in Images]
    x_vals = np.arange(len(data))
    # (1) show the 5-95 by the notches, *dont* show any outliers (too slow)
    # (2) show the median line 
    box = plt.boxplot(data,positions=x_vals,whis=[5,95],sym='',
                      patch_artist=True)
    plt.xticks(x_vals, labels,rotation=0,)
    PlotUtilities.lazyLabel("","Height (nm)",
                            "Comparison of Image Distributions")
    # set the alpha and colors
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.3)
    # finally, plot the mean and standard deviations as a functo
    PlotUtilities.savefig(fig,OutBase+ "Grid_" + Name+".pdf")
            

def InitialDebugging():
    """
    Makes the plots for the early 2016-10 images of surfaces and tips
    """
    surface_kwargs = dict(vmin=-1, vmax=4)
    tip_kwargs = dict(vmin=-10,vmax=18)
    Base = "/Volumes/group/4Patrick/Reports/" + \
           "2016-10-20-tip-and-surface-height-distributions/"
    files_tips =[
        ["1_EtchedLong.pxp","Etched Long",tip_kwargs],
        ["9-221-2016-functiuonalized-9-19-tStrept-in-pbs-ph7.4-with-mini.pxp","'Good' T-Strept Long",tip_kwargs],
        ["2016-10-14-mini-pbs-t-strept-tip-batch-10-10-2016_1%_attachment.pxp","'Bad'(?) T-Strept Long",tip_kwargs] ]
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

def SilanePegMalemideDebugging():
    """
    Makes the plots for the Silane-Peg-Maleimide and UVChamber/Lamp comparisons
    from 10-25
    """
    common_kwargs = dict(vmin=-1, vmax=6)
    Base = "/Volumes/group/4Patrick/Reports/" + \
           "2016-10-25-spm-comparison/"
    files_tips =[
         ["2016-10-18-mini-pbs-PEG-azide-surface-batch-1x_peg_0p15mg_mL_10-17-2016_0%_attachment_with_proteins_possibly-tips-fault.pxp",
          "SPA on Glass, UV Chamber",common_kwargs],
        ["2016-10-25-spm-on-kohd-glass-in-pbs-imaged-by-mini.pxp",
         "SPM on Glass, UV Chamber",common_kwargs],
        ["2016-10-25-spm-on-biol-uv-chamber-not-lamp-in-pbs-imaged-by-mini.pxp",
         "SPM Bio-L, UV Chamber",common_kwargs],
        ["2016-10-25-spm-on-biol-uv-lamp-not-chamber-in-pbs-imaged-by-mini.pxp",
         "SPM Bio-L, UV Lamp",common_kwargs]
          ]
    MakeGridPlot(files_tips,[-5,40],Base,"SurfaceGrid",figsize=(24,16),
                 ImageFunc=lambda m_list: m_list[-1])

def HigherMwPEG():
    """
    Makes the plots for the 3kDA Silane-Peg-Maleimide from 10/26
    """
    common_kwargs = dict(vmin=-1, vmax=10)
    Base = "/Volumes/group/4Patrick/Reports/" + \
           "2016-10-26-3kDa/"
    files_tips =[
         ["2016-10-18-mini-pbs-PEG-azide-surface-batch-1x_peg_0p15mg_mL_10-17-2016_0%_attachment_with_proteins_possibly-tips-fault.pxp",
          "SPA on Glass",common_kwargs],
        ["2016-10-25-spm-on-kohd-glass-in-pbs-imaged-by-mini.pxp",
         "600Da SPM on Glass",common_kwargs],
        ["2016-10-26-3kDa-spm-on-kohd-glass-in-pbs-imaged-by-mini.pxp",
         "3kDa SPM on Glass",common_kwargs],
        ["2016-10-25-spm-on-biol-uv-chamber-not-lamp-in-pbs-imaged-by-mini.pxp",
         "600Da SPM on Bio-L",common_kwargs],
        ["2016-10-26-3kDa-spm-on-bio-long-in-pbs-imaged-by-mini.pxp",
         "3kDa SPM on Bio-L",common_kwargs],
        ["2016-10-28-imaging-spm-3kDa-bio-l-low-attachments-on-DNA-surface-in-pbs-or-dried-dna-surface.pxp",
         "3kDa SPM-Strept on Bio-L",common_kwargs]
          ]
    MakeGridPlot(files_tips,[-5,40],Base,"SurfaceGrid",figsize=(28,22),
                 ImageFunc=lambda m_list: m_list[-1],
                 hist_kwargs=dict())

    
def CleaningConditions():
    """
    Makes the plots for the cleaning conditions trials from 10/28
    """
    common_kwargs = dict(vmin=-1, vmax=10)
    Base = "/Volumes/group/4Patrick/Reports/2016-10-28-cleaning-conditions/"
    files_tips =[
         ["2016-10-28-spm-600Da-on-glass-pos-ctrl.pxp",
          "600 SPM, Glass, Double UV",common_kwargs],
        ["2016-10-25-spm-on-biol-uv-chamber-not-lamp-in-pbs-imaged-by-mini.pxp",
         "600 SPM, Bio-L, 'Usual' Clean",common_kwargs],
        ["2016-10-28-spm-600Da-on-uv-ozoned-biolong-cleaned-before-and-after-etch-with-solvent-rinse-after-etch.pxp","SPM, Bio-L, Double UV + Solvent",
         common_kwargs],
        ["2016-10-28-spm-600Da-on-uv-ozoned-biolong-cleaned-before-and-after-etch-no-solvent-after.pxp","600 SPM, Bio-L, Double UV",common_kwargs]
        ]
    MakeGridPlot(files_tips,[-5,40],Base,"SurfaceGrid",figsize=(20,16),
                 ImageFunc=lambda m_list: m_list[-1],
                 hist_kwargs=dict())
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    CleaningConditions()
    HigherMwPEG()
    SilanePegMalemideDebugging()
    InitialDebugging()
    
if __name__ == "__main__":
    run()
