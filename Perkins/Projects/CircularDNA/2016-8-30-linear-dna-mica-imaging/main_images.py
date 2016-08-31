# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from skimage.io import imread
import matplotlib.cm as cm
from GeneralUtil.python import PlotUtilities
## sci-kit image tools for image segmentation. See:
# scikit-image.org/docs/dev/user_guide/tutorial_segmentation.html
# watershed provides segmentation
from skimage.morphology import watershed
# sobel alborithm computed amplitude of gradient, good for elevation mapping
from skimage.filters import sobel
# used for labelling and hole filling
from skimage.measure import label,regionprops
from skimage.morphology import skeletonize,medial_axis
import matplotlib.patches as mpatches

import cv2

def RegionArcLength(RegionProperties):
    coords = RegionProperties.coords
    return cv2.arcLength(coords,closed=False)

def run():
    """
    Simple demo for reading in and segmenting an image
    """
    ImagePath = "./Image0007HtR_Small.tif"
    # read in the image and background correct it
    Image = imread(ImagePath,as_grey=True)
    bg = np.median(Image)
    CorrectedImage = Image - bg
    # get the threshholding markers for this image. We are very generous with
    # what we are sure is background and restrictive on what
    # is considered signal
    qLow,qHigh = np.percentile(CorrectedImage,[5,95])
    LowIntensity = qLow
    HighIntensity = qHigh
    # determine 'markers' for the data and background
    markers = np.zeros_like(CorrectedImage)
    markers[CorrectedImage <= LowIntensity] = 1
    markers[CorrectedImage >= HighIntensity] = 2
    # create an elevation map to determine edges
    elevation_map = sobel(CorrectedImage)
    # segment based on the markers and elevation map
    segmentation = watershed(elevation_map, markers)
    # skeletonize by going to 0/1 representaiton (currently 1/2)
    skeleton,dist = medial_axis(segmentation - 1,return_distance=True)
    print(dist)
    # get the labelled regions (this consists of contiguous pieces of DNA)
    labels = label(skeleton)
    # get the properties of the labells (this gives us the chain statistics)
    chain_properties = regionprops(labels)
    # get the arc lengths (contour lengths) of the regions
    arc_lengths = np.array([RegionArcLength(r) for r in  chain_properties])
    PixelsPerNanometer = 2e3/1024
    print("\n".join("{:.1f}nm".format(l)
                    for l in arc_lengths*PixelsPerNanometer))
    fig = PlotUtilities.figure(dpi=600)
    # get the arc length of each chain
    NumRows = 5
    plt.subplot(NumRows,2,1)
    plt.imshow(Image, cmap=cm.gray)
    plt.subplot(NumRows,2,3)
    plt.imshow(CorrectedImage, cmap=cm.gray)
    plt.subplot(NumRows,2,4)
    plt.hist(CorrectedImage.ravel(),bins=100)
    plt.axvline(LowIntensity)
    plt.axvline(HighIntensity)
    plt.subplot(NumRows,2,5)
    plt.imshow(elevation_map,cmap=cm.gray)
    plt.subplot(NumRows,2,7)
    plt.imshow(segmentation,cmap=cm.gray)
    # draw and indicate the skeletons
    ax = plt.subplot(NumRows,2,8)
    plt.imshow(skeleton,cmap=cm.Greys)
    colors = ['r','g','b','k','m','y']
    for i,region in enumerate(chain_properties):
        # draw rectangle around region coins
        minr, minc, maxr, maxc = region.bbox
        color = colors[i % len(colors)]
        rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
                                  fill=False, edgecolor=color, linewidth=1)
        ax.add_patch(rect)
    PlotUtilities.savefig(fig,"image.png")

if __name__ == "__main__":
    run()
