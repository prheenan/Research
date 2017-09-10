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
sys.path.append("../../../../../../../")
from Research.Perkins.AnalysisUtil.Images import  ImageUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from GeneralUtil.python import GenUtilities,PlotUtilities,CheckpointUtilities
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import gaussian
from skimage.morphology import skeletonize
from skimage import measure

class transform:
    def __init__(self,name,function,imshow_kw=dict(cmap=plt.cm.afmhot)):
        self.name = name
        self.function = function
        self.imshow_kw = imshow_kw

def subtract_background(image,deg=2,**kwargs):
    image = image.height
    to_ret = image.copy().T
    shape = image.shape
    coords = np.arange(shape[1])
    coeffs = np.array(np.polyfit(x=coords,y=image,deg=deg,**kwargs))
    n_rows = shape[0]
    for i in range(n_rows):
        to_ret[i,:] -= np.polyval(coeffs[:,i],x=coords)
    return to_ret.T

def realistic_min_max(image,q_low=90,q_high=99.9):
    data = image.height_nm()
    low,high = np.percentile(data.flatten(),[q_low,q_high])
    return low,high
    
def plot_with_background_corrected(args,imshow_kw_list=None):
    n = len(args)
    fig = PlotUtilities.figure((3.5,2*n))
    for i,a in enumerate(args):
        ax = plt.subplot(n,1,(i+1))
        vmin,vmax = realistic_min_max(a)
        if imshow_kw_list is not None:
            m_list = imshow_kw_list[i]
        else:
            m_list = dict()
        imshow_kwargs = dict(vmin=vmin,vmax=vmax,**m_list)
        im = ImageUtil.make_image_plot(a,pct=50,imshow_kwargs=imshow_kwargs)
        if (i == 0):
            ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig)
        if (i < n-1):
            PlotUtilities.xlabel("",ax=ax)
            PlotUtilities.no_x_label(ax=ax)
            ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig,add_space_only=True)
        else:
            ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig,add_space_only=True)
    return fig

def _safe_apply(images,f):
    to_ret = []
    for ex in images:
        tmp = copy.deepcopy(ex)
        ret = f(ex)
        assert type(ret) is np.ndarray , "{:s} didn't return array".format(f)
        tmp.height = ret
        to_ret.append(tmp) 
    return to_ret 



def threshold(im,threshold_nm,rel_pct=50):
    """
    Returns: im, height zeroed where height-rel_pct_of_height < threshold_nm
    """
    height_rel = im.height_nm() 
    height_rel -= np.percentile(height_rel,rel_pct)
    zero_idx = np.where(height_rel < threshold_nm)
    height_new = copy.deepcopy(im.height)
    height_new[zero_idx] = 0
    return height_new

def binarize(image):
    """
    binarizes a single image: set to 1 where the image is non-zero
    """
    binary = copy.deepcopy(image)
    binary[binary > 0] = 1
    return binary

def correct_background(images,**kw):
    """
    See: threshold_images, except subtracts the AFM background 
    """
    return _safe_apply(images,lambda x: subtract_background(x,**kw))



def blur_images(images,sigma=1,**kw):
    """
    See: thresholdimages, except adds a gaussian blur with sigma
    """
    return _safe_apply(images,lambda x: gaussian(x.height,sigma=sigma,**kw))

def threshold_images(images,threshold_nm=0.2):
    """
    thresholds the heights for each of images. pass by copy
    """
    return _safe_apply(images,lambda x: threshold(x,threshold_nm))

def binarize_images(images):
    """
    binarizes the heights for each of images. pass by copy
    """
    return _safe_apply(images,lambda x: binarize(x.height))

def skeletonize_images(images):
    """
    returns: the skeletonized version of the (assumed already binary) images
    """
    return _safe_apply(images,lambda x: skeletonize(x.height))

def label_images(images):
    """
    returns: the labelled versions (ie: connected components) of the (assumed
             skeletonized) images
    """
    return _safe_apply(images,lambda x: measure.label(x.height,background=0))

def cache_images(cache_dir,func,**kw):
    """
    either caches or re-func to get an image transformaiton

    Args:
        cache_dir: where the cache is, or will be created
        func: functor (no arguments), re-reads everything if needed
        **kw: passed to CheckpointUtilities.multi_load
    Returns:
        cached imags
    """
    return CheckpointUtilities.multi_load(cache_dir,load_func=func,
                                          name_func=FEC_Util.name_func,**kw)



    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "./data/"
    cache_dir_fmt = "./cache_{:s}/"
    cache_dir_raw = cache_dir_fmt.format("raw")
    out_dir = in_dir
    force_def = dict(force=True)
    GenUtilities.ensureDirExists(out_dir)
    images = ImageUtil.cache_images_in_directory(pxp_dir=in_dir,
                                                 cache_dir=cache_dir_raw,
                                                 limit=2)
    corrected_dir = cache_dir_fmt.format("corrected")
    images = [images[1]]
    transforms = [transform("corrected",correct_background),
                  transform("gaussian",blur_images),
                  transform("threshold",threshold_images),
                  transform("binarize",binarize_images),
                  transform("skeletonize",skeletonize_images),
                  transform("label",label_images,dict(cmap=plt.cm.spectral))]
    last = images
    all_transforms = [last]
    for tx in transforms:
        tmp_dir = cache_dir_fmt.format(tx.name)
        last = cache_images(tmp_dir,func = lambda: tx.function(last),
                            **force_def)
        all_transforms.append(last)
    last_dir = tmp_dir
    # subtract the linear backround from each, save to a new cache 
    for i in range(len(images)):
        pipeline = [x[i] for x in all_transforms]
        # add in the keywords for the first image...
        imshow_kw = [dict(cmap=plt.cm.afmhot)] + \
                    [tx.imshow_kw for tx in transforms]
        fig = plot_with_background_corrected(pipeline,imshow_kw_list=imshow_kw)
        img_name = FEC_Util.name_func(i,images[i])
        out_name = "{:s}{:s}.png".format(last_dir,img_name)
        PlotUtilities.savefig(fig,out_name)

    
    
                                                 
if __name__ == "__main__":
    run()
