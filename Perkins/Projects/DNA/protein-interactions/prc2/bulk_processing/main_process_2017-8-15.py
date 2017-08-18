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
sys.path.append("../../../../../../../")
from Research.Perkins.AnalysisUtil.Images import  ImageUtil
from GeneralUtil.python import GenUtilities,PlotUtilities
from mpl_toolkits.axes_grid1 import make_axes_locatable

def make_image_plot(im,imshow_kwargs,pct=1):
    # offset the data
    im_height = im.height_nm()
    min_offset = np.median(im_height)
    im_height -= min_offset
    range_microns = im.range_meters * 1e6
    to_ret = plt.imshow(im_height.T,extent=[0,range_microns,0,range_microns],
                        interpolation='bicubic',**imshow_kwargs)
    PlotUtilities.tom_ticks()
    micron_str = PlotUtilities.upright_mu("m")
    PlotUtilities.lazyLabel(micron_str,micron_str,"",
                            tick_kwargs=dict(direction='out'))    
    return to_ret                             

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "./data/"
    cache_dir = "./cache/"
    out_dir = in_dir
    GenUtilities.ensureDirExists(out_dir)
    images = ImageUtil.cache_images_in_directory(pxp_dir=in_dir,
                                                 cache_dir=cache_dir)
    # note: vmin/vmax are in nm (as is height)                                                  
    vmin_nm = 0.0
    vmax_nm = vmin_nm + 1.5
    imshow_kwargs = dict(vmin=vmin_nm,vmax=vmax_nm,cmap = plt.cm.Greys_r)    
    for im in images:
        src_file = im.SourceFilename()
        # path to save the image out 
        full_out_path = "{:s}{:s}-{:s}.png".format(out_dir,src_file,im.Name())        
        fig = PlotUtilities.figure((3.5,4))
        ax = plt.subplot(1,1,1)
        im = make_image_plot(im,imshow_kwargs,pct=1) 
        # make a separate axis for the colorbar 
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        PlotUtilities.colorbar("Height (nm)",fig=fig,
                               bar_kwargs=dict(mappable=im,cax=cax))
        PlotUtilities.savefig(fig,full_out_path,bbox_inches='tight')

if __name__ == "__main__":
    run()
