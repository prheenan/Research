# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,matplotlib

sys.path.append("../../../../../../../../")

from Research.Perkins.AnalysisUtil.Images import ImageUtil
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import GenUtilities,PlotUtilities
## sci-kit image tools for image segmentation. See:
# scikit-image.org/docs/dev/user_guide/tutorial_segmentation.html
# watershed provides segmentation
from skimage.morphology import watershed
from skimage.segmentation import active_contour
# sobel alborithm computed amplitude of gradient, good for elevation mapping
from skimage.filters import sobel
# used for labelling and hole filling
from skimage.measure import label,regionprops
from skimage.morphology import skeletonize,medial_axis
import matplotlib.patches as mpatches
from skimage.morphology import medial_axis, skeletonize, skeletonize_3d

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base_dir = "./data/"
    out_dir = "./out/"
    cache_dir = "./cache/"
    GenUtilities.ensureDirExists(cache_dir)
    GenUtilities.ensureDirExists(out_dir)    
    images = ImageUtil.cache_images_in_directory(base_dir,cache_dir)
    # only get images which are two microns
    images = [i for i in images if abs(i.range_meters - 2e-6) < 1e-9]
    for i in images:
        # manually set the images size so that it is 'all image' 
        dpi = 600
        figsize = np.array(i.height.shape) / dpi
        fig = plt.figure(frameon=False,dpi=dpi)
        fig.set_size_inches(figsize)
        ax = plt.Axes(fig,[0,0,1,1])
        ax.set_axis_off()
        fig.add_axes(ax)
        height_nm_rel = i.height_nm_rel()
        vmin,vmax = np.percentile(height_nm_rel,[80,99.5])
        ax.imshow(height_nm_rel,cmap=plt.cm.Greys_r,vmin=vmin,vmax=vmax,
                  aspect='normal')
        size_microns = i.range_meters * 1e6                
        size_str = "_{:.2g}microns".format(size_microns)
        save_name =out_dir + i.SourceFilename() + i.Name() + size_str + ".tiff"
        fig.savefig(save_name,dpi=dpi)


if __name__ == "__main__":
    run()
