# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
from Research.Perkins.AnalysisUtil.Images import ImageUtil

import matplotlib.gridspec as gridspec
from scipy.stats import norm

def ReadImageAsObject(file_name):
    """
    Convenience wrapper; reads the file as an image
    
    Args:
        file_name: See ReadImage
    Returns:
        List of SurfaceImage objects present in the file
    """
    images = PxpLoader.ReadImage(file_name)
    return [PxpLoader.SurfaceImage(ex) for ex in images]
def MakePlot(SurfaceImage,label,OutPath,**kwargs):
    """
    Makes a simple plot of the desired distribution
    """
    # make plots
    fig = PlotUtilities.figure(figsize=(10/1.5,16/1.5))
    ax = plt.subplot(2,1,1)
    ImageUtil.PlotImage(SurfaceImage,**kwargs)
    PlotUtilities.lazyLabel(r"Microns",r"Microns",label)
    PlotUtilities.colorbar("Height (nm)")
    ax = plt.subplot(2,1,2)
    bins_tmp = ImageUtil.PlotImageDistribution(SurfaceImage)
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
        ImageUtil.PlotImage(im,**kwargs)
        if (i ==0 ):
            PlotUtilities.lazyLabel(r"Microns",r"Microns",label)
            PlotUtilities.colorbar("Height (nm)")
        else:
            PlotUtilities.lazyLabel(r"","",label)
            PlotUtilities.colorbar("")
        ax = plt.subplot(gs[1,i])
        bins_tmp = ImageUtil.PlotImageDistribution(im,color=colors[color_idx])
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
        ImageUtil.PlotImageDistribution(im,label=labels[i],
                                        color=colors[color_idx],
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
        ["9-221-2016-functiuonalized-9-19-tStrept-in-pbs-ph7.4-with-mini.pxp",
         "'Good' T-Strept Long",tip_kwargs],
        ["2016-10-14-mini-pbs-t-strept-tip-batch-10-10-2016_1%_attachment.pxp",
         "'Bad'(?) T-Strept Long",tip_kwargs] ]
    files_surfaces = [
        ["2016-10-6-mini-pbs-koh-clened-glass-batch-10-5-2016.pxp","KOH",
         surface_kwargs],
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

def PositiveControls():
    """
    Makes the plots for the positive controls (bio-L, t-strept and peg) we have
    from week of 11/7/2016
    """
    common_kwargs = dict(vmin=-1, vmax=10)
    Base = "/Volumes/group/4Patrick/Reports/2016-11-14-working-positive-control-tstrept-bio-l/"
    files_tips =[
         ["2016-11-8-spm-3kDa-11-7-batch-biol-biol-mini-imaged.pxp",
          "3kDa PEG, no t-strept",common_kwargs],
        ["2016-11-8-tstrep-tip-decent-attachment-made-11-1-2016-biol-mini-imaged.pxp", "3kDa, t-strept, 1 week old",common_kwargs],
        ["2016-11-8-tstrept-004-good-attachment-11-7-batch-biol-biol-mini-imaged.pxp","3kDa, t-strept, new, 1",
         common_kwargs],
        ["2016-11-9-tstrept-005-ty-multiple-freeze-thaw-good-attachment-11-7-batch-biol-biol-mini-imaged.pxp","3kDa, t-stept, new, 2",common_kwargs],
        ["2016-11-9-tstrept-005-patrick-single-freeze-thaw-good-attachment-11-7-batch-biol-biol-mini-imaged.pxp","3kDa, t-strept, new, 3",common_kwargs]]
    MakeGridPlot(files_tips,[-5,40],Base,"SurfaceGridPosCtrl",figsize=(25,16),
                 ImageFunc=lambda m_list: m_list[-1],
                 hist_kwargs=dict())

def ImageFunction(m_list):
    assert len(m_list) > 0, "Didn't find any images..."
    return m_list[-1]

def KOH_Testing():
    common_kwargs = dict(vmin=-2, vmax=3)
    Base = "/Volumes/group/4Patrick/Reports/2017-2-17-koh-comparisons/"
    files_tips =[
         ["2017-1-27-stephen-2.5M-koh-glass-in-pbs.pxp",
          "KOH, 95% ethanol",common_kwargs],
        ["2017-1-27-traditional-5.7M-koh-glass-in-pbs.pxp",
         "KOH, 100% ethanol",common_kwargs]]
    MakeGridPlot(files_tips,[-5,40],Base,"SurfaceGridPosCtrl",figsize=(25,16),
                 ImageFunc=ImageFunction,
                 hist_kwargs=dict(),Force=False)
    
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    KOH_Testing()
    PositiveControls()
    CleaningConditions()
    HigherMwPEG()
    SilanePegMalemideDebugging()
    InitialDebugging()
    
if __name__ == "__main__":
    run()
