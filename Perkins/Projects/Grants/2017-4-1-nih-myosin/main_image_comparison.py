# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import PlotUtilities,CheckpointUtilities,GenUtilities
from Research.Perkins.AnalysisUtil.Images import ImageUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

import matplotlib.gridspec as gridspec
from scipy.stats import norm
import matplotlib
from matplotlib_scalebar.scalebar import ScaleBar


    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = FEC_Util.default_data_root()
    input_directory = base + \
        "4Patrick/CuratedData/Grants/2017-4-1-nih-myosin/igor/cropped"
    files = GenUtilities.getAllFiles(input_directory,ext=".ibw")
    # read all the ibw files in and cache them
    images = []
    func = PxpLoader.read_ibw_as_image
    for i,f in enumerate(files):
        cache_file = "./{:d}.pkl".format(i)
        tmp = CheckpointUtilities.getCheckpoint(cache_file,func,False,f)
        images.append(tmp)
    pcts = [np.percentile(i.height_nm_rel(),[75,99]) for i in images]
    vmin = np.min(pcts)
    vmax = np.max(pcts)
    n_plots = len(images)
    fig,axes = plt.subplots(figsize=(4*n_plots,3),nrows=1,ncols=n_plots)        
    ax = []
    for i,image in enumerate(images):  
        ax = axes.flat[i]
        range_plot_nanometers = 1000 * image.range_microns()   
        vmin_dict = dict(vmin=vmin,vmax=vmax)
        im = ImageUtil.PlotImage(image,cmap=plt.cm.afmhot,fix_extent=False,
                                 ax=ax,range_plot=range_plot_nanometers,
                                 **vmin_dict)
        PlotUtilities.FormatImageAxis()
    fig.colorbar(im,ax=list(axes),label="Height (nm)")
    # for some reason, need to reformat
    for ax,image in zip(axes,images):
        PlotUtilities.FormatImageAxis(ax=ax)   
        pixel_size_meters = image.pixel_size_meters
        pixel_size_nanometers = pixel_size_meters * 1e9
        scalebar = ScaleBar(pixel_size_nanometers,'nm',box_alpha=0.7)
        ax.add_artist(scalebar)     
    fudge_x = 0.1
    fudge_y = 0.02
    subplots_adjust=dict(top=1-fudge_y,bottom=fudge_y,
                         left=fudge_y,right=1-fudge_y,
                         hspace=fudge_y,wspace=fudge_x)    
    PlotUtilities.savefig(fig,"{:d}.png".format(i),tight=False,
                          subplots_adjust=subplots_adjust)
    

        
if __name__ == "__main__":
    run()
