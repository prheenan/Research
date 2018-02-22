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
import matplotlib.colors as colors

def _image_helper(im,fig,ax,imshow_kwargs=None,colorbar_kwargs=dict(),
                  vmin_nm=-0.1,vmax_nm=4,**kw):
    """
    Utility function for creating reasonable images
    
    Args:
        im: the PxpLoader.SurfaceImage object   
        fig,ax: the figure and axes  to use (important for colorbar)
        imshow_kwargs: passed to ImageUtil.make_image_plot
        colorbar_kwargs: passed to smart_colorbar
    """
    if (imshow_kwargs is None):                                
        imshow_kwargs = dict(vmin=vmin_nm,vmax=vmax_nm,cmap = plt.cm.afmhot,
                             norm=None)
    try:
        im.height = ImageUtil.subtract_background(im)
    except ValueError as e:
        print(e)
    im = ImageUtil.make_image_plot(im,imshow_kwargs,pct=50,**kw)               
    ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig,**colorbar_kwargs)     


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    pass

if __name__ == "__main__":
    run()