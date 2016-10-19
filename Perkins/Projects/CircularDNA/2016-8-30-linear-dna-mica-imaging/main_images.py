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
from skimage.segmentation import active_contour
# sobel alborithm computed amplitude of gradient, good for elevation mapping
from skimage.filters import sobel
# used for labelling and hole filling
from skimage.measure import label,regionprops
from skimage.morphology import skeletonize,medial_axis
import matplotlib.patches as mpatches
from skimage.filters import threshold_adaptive
# cv2 for arc length
from cv2 import arcLength

# XXX to check out:
# fitting snakes to data:
#scikit-image.org/docs/dev/auto_examples/edges/plot_active_contours.html

def RegionArcLength(RegionProperties):
    coords = RegionProperties.coords
    return arcLength(coords,closed=False)

def run():
    """
    Simple demo for reading in and segmenting an image
    """
    ImagePath = "./Image0007HtR_Small.tif"
    PixelsPerNanometer = 2e3/1024
    MinFeatureSizeNanometer = 10
    MinPixels = int(np.ceil(MinFeatureSizeNanometer/PixelsPerNanometer))
    MinPixels= MinPixels + 1 if (MinPixels % 2 == 0) else MinPixels
    # read in the image and background correct it
    Image = imread(ImagePath,as_grey=True)
    # denoising, see
    # scikit-image.org/docs/dev/auto_examples/filters/plot_denoise.html
    qLow,qHigh = np.percentile(Image,[90,99])
    LowIntensity = qLow
    HighIntensity = qHigh
    CorrectedImage =  Image.copy()
    CorrectedImage[CorrectedImage < LowIntensity] = 0
    CorrectedImage = CorrectedImage/np.max(CorrectedImage)
    # get the threshholding markers for this image. We are very generous with
    # what we are sure is background and restrictive on what
    # is considered signal
    # determine 'markers' for the data and background
    markers = np.zeros_like(CorrectedImage)
    markers[CorrectedImage <= LowIntensity] = 1
    markers[CorrectedImage >= HighIntensity] = 2
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
    print("\n".join("{:.1f}nm".format(l)
                    for l in arc_lengths*PixelsPerNanometer))
    fig = PlotUtilities.figure(dpi=600)
    # get the arc length of each chain
    NumRows = 5
    ImageFig = lambda : PlotUtilities.figure()
    IncrementName = lambda fig,x: \
                    PlotUtilities.savefig(plt.gcf(),"Image{:d}.png".format(x))
    ImageFig()
    plt.imshow(Image, cmap=cm.gray)
    IncrementName(fig,0)
    ImageFig()
    plt.imshow(CorrectedImage, cmap=cm.gray)
    IncrementName(fig,1)
    ImageFig()
    plt.imshow(binary,cmap=cm.gray)
    IncrementName(fig,2)
    ImageFig()
    plt.imshow(segmentation,cmap=cm.gray)
    IncrementName(fig,3)
    # draw and indicate the skeletons
    ImageFig()
    ax = plt.subplot(1,1,1)
    plt.imshow(skeleton,cmap=cm.Greys)
    colors = ['r','g','b','k','m','y']
    for i,region in enumerate(chain_properties):
        # draw rectangle around region coins
        minr, minc, maxr, maxc = region.bbox
        color = colors[i % len(colors)]
        rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
                                  fill=False, edgecolor=color, linewidth=1)
        ax.add_patch(rect)
    IncrementName(fig,4)

if __name__ == "__main__":
    run()
