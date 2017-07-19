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

sys.path.append("../../../../../../../")

from IgorUtil.PythonAdapter import PxpLoader
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
from skimage.filters import threshold_adaptive
import cv2
from skimage.filters import gaussian
from skimage.morphology import medial_axis, skeletonize, skeletonize_3d


def RegionArcLength(RegionProperties):
    coords = RegionProperties.coords
    return cv2.arcLength(coords,closed=False)

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    example_file = "./prh_tmp_im_28.ibw"
    surface_image = PxpLoader.read_ibw_as_image(example_file)
    height_nm = surface_image.height_nm()
    min_height = np.min(height_nm)
    img = ((height_nm - min_height)/(np.max(height_nm)-min_height))
    img = img[150:250,375:500]
    img = gaussian(img,1)
    # write down a line connecting the known endpoints
    x = list(np.linspace(32,49,num=50)) + \
        list(np.linspace(49,49,num=50)) + \
        list(np.linspace(49,117,num=50))
    y = list(np.linspace(67,81,num=50)) + \
        list(np.linspace(81,55,num=50)) + \
        list(np.linspace(55,20,num=50))

    init = np.array([x, y]).T
    #snake = active_contour(img, init, bc='fixed',max_px_move=2,
    #                       alpha=0.1, beta=1.0, w_line=5, w_edge=-10, gamma=1,
     #                      convergence=1e-9,max_iterations=int(1e4))
    # Compute the medial axis (skeleton) and the distance transform
    skel, distance = medial_axis(img, return_distance=True)
    dist_on_skel = distance * skel
    plt.imshow(img, cmap=plt.cm.gray, interpolation='nearest')
    plt.imshow(dist_on_skel, cmap=plt.cm.viridis, interpolation='nearest')
    plt.show()
    exit(1)
    # Distance to the background for pixels of the skeleton

    # denoising, see
    # scikit-image.org/docs/dev/auto_examples/filters/plot_denoise.html
    qLow,qHigh = np.percentile(img,[60,95])
    LowIntensity = qLow
    HighIntensity = qHigh
    CorrectedImage =  img.copy()
    CorrectedImage[CorrectedImage < LowIntensity] = 0
    CorrectedImage[CorrectedImage >= HighIntensity] = 1
    spline_x = snake[:, 0]
    spline_y = snake[:, 1]
    diff_spline_x_sq = np.diff(spline_x)**2
    diff_spline_y_sq = np.diff(spline_y)**2
    diff_spline_sq = np.diff(spline_x)**2 + np.diff(spline_y)**2
    diff_spline = np.sqrt(diff_spline_sq)
    arc_length_pixels = sum(diff_spline)
    arc_length_meters = arc_length_pixels * surface_image.pixel_size_meters
    arc_length_microns = arc_length_meters * 1e6
    print("Length of snake is {:.3g}um".format(arc_length_microns))
    plt.imshow(CorrectedImage,cmap=plt.cm.gray)
    plt.plot(init[:, 0], init[:, 1], '--r', lw=3)
    plt.plot(spline_x,spline_y, '-b', lw=3)
    plt.show()
    
    exit(1)


    # get the threshholding markers for this image.
    markers = np.zeros(CorrectedImage.shape)
    markers[CorrectedImage <= LowIntensity] = 1
    markers[CorrectedImage > LowIntensity] = 2
    # create an elevation map to determine edges
    elevation_map = sobel(CorrectedImage)
    # segment based on the markers and elevation map
    segmentation = watershed(elevation_map, markers)
    binary =  segmentation - 1
    # skeletonize by going to 0/1 representaiton (currently 1/2)
    skeleton = skeletonize(binary)
    # get the labelled regions (this consists of contiguous pieces of DNA)
    labels = label(skeleton)
    # get the properties of the labells (this gives us the chain statistics)
    chain_properties = regionprops(labels)
    # get the arc lengths (contour lengths) of the regions
    arc_lengths = np.array([RegionArcLength(r) for r in  chain_properties])
    plt.imshow(skeleton)
    for i,c in enumerate(chain_properties):
        if (arc_lengths[i] > 50):
            plt.plot(*c.centroid,color='r',marker='o')
    plt.show()


if __name__ == "__main__":
    run()
