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
from matplotlib.ticker import AutoMinorLocator


    
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
    # rotate the first image
    example = images[0]
    example.rotate(angle_degrees=140,reshape=False)
    pcts = [np.percentile(i.height_nm_rel(),[60,100]) for i in images]
    vmin = np.min(pcts)
    vmax = np.max(pcts)
    n_plots = len(images)
    width = 3.5*n_plots
    fig,axes = plt.subplots(figsize=(width,3),nrows=1,ncols=n_plots)        
    ax = []
    font_common=dict(fontname="Arial")
    for i,image in enumerate(images):  
        ax = axes.flat[i]
        range_plot_nanometers = 1000 * image.range_microns()   
        vmin_dict = dict(vmin=vmin,vmax=vmax)
        im = ImageUtil.PlotImage(image,cmap=plt.cm.afmhot,fix_extent=False,
                                 ax=ax,range_plot=range_plot_nanometers,
                                 **vmin_dict)
    cax = fig.add_axes([0.91,0.1,0.03,0.80])
    cbar = fig.colorbar(im,ax=list(axes),label="Height (nm)",cax=cax,
                        ticks=[0,1,2,3])
    cbar.ax.minorticks_on()              
    cbar.ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    common_tick = dict(direction='in')
    ax.tick_params(which='minor',length=8,width=2,**common_tick)
    ax.tick_params(which='major',length=15,width=4,**common_tick)
    # for some reason, need to reformat
    for i,(ax,image )in enumerate(zip(axes,images)):
        PlotUtilities.FormatImageAxis(ax=ax)   
        pixel_size_meters = image.pixel_size_meters
        pixel_size_nanometers = pixel_size_meters * 1e9
        scalebar = ScaleBar(pixel_size_nanometers,'nm',box_alpha=0,
                            location=(1),color='w',length_fraction=0.3)
        ax.add_artist(scalebar)   
    fudge_x = 0.01
    fudge_y = 0.005
    subplots_adjust=dict(top=1-fudge_y,bottom=fudge_y,
                         left=fudge_y,right=0.9,
                         hspace=fudge_y,wspace=fudge_x)    
    axis_func = lambda x: x[:-1]                             
    PlotUtilities.label_tom(fig,loc=(0.05,0.93),color='w',fontsize=15,
                            axis_func=axis_func,**font_common)                        
    PlotUtilities.savefig(fig,"{:d}.png".format(i),tight=True,
                          subplots_adjust=subplots_adjust)    

        
if __name__ == "__main__":
    run()
