# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,copy

sys.path.append("../")
sys.path.append("../../../../../../../../../")
from Research.Perkins.AnalysisUtil.Images import ImageUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

from GeneralUtil.python import CheckpointUtilities,PlotUtilities
from Util import Processing
from skimage.morphology import skeletonize
from skimage.filters import gaussian
from skimage.measure import label,regionprops

def skel(image,thresh):  
    binary = np.zeros_like(image)
    cond_ge = (image >= thresh)
    binary[np.where(cond_ge)] = 1
    skeleton = skeletonize(binary)
    return skeleton

def filter(images,threshold_m,blur_sigma_px):
    for i_original in images:
        i = copy.deepcopy(i_original)
        height_m = i.height
        height_m_rel = height_m - np.median(height_m)
        height_m_rel_corr = ImageUtil._subtract_array_background(height_m_rel)
        # blur the image, to get rid of pixel noise 
        height_m_blurred = gaussian(height_m_rel_corr,sigma=blur_sigma_px)
        skeleton = skel(height_m_blurred,threshold_m)
        label_image = label(skeleton)
        props = regionprops(label_image=label_image,cache=True,
                             intensity_image=height_m_rel_corr)
        to_ret = Processing.TraceInfo(raw_image=i,trace_properties=props)
        yield to_ret
        
def run():
    base_dir = Processing.cache_select_images()
    threshold = 0.3e-9
    blur_sigma_px = 1
    cache_dir = Processing.cache_traces()
    images = CheckpointUtilities.lazy_multi_load(base_dir)
    load_func = lambda: filter(images,threshold,blur_sigma_px)
    e = CheckpointUtilities.multi_load(cache_dir=cache_dir,load_func=load_func,
                                       name_func=FEC_Util.fec_name_func)
    for tmp in e:
        height_rel = tmp.raw_image.height
        height_rel -= np.median(height_rel)
        fig = PlotUtilities.figure()        
        ax1 = plt.subplot(2,1,1)
        plt.imshow(height_rel,vmax=2e-9,vmin=0)
        PlotUtilities.FormatImageAxis(ax1)
        ax2 = plt.subplot(2,1,2)
        plt.imshow(tmp.label_image,cmap=plt.cm.spectral,interpolation='nearest')
        PlotUtilities.FormatImageAxis(ax2)
        out_name = cache_dir + FEC_Util._fec_name_func(tmp)+".jpeg"
        PlotUtilities.savefig(fig,out_name)
        
if __name__ == "__main__":
    run()

    