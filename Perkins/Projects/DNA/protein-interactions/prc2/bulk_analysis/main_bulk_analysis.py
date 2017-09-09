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

def subtract_background(image,deg=2,**kwargs):
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
    return dict(vmin=low,vmax=high)
    
def plot_with_background_corrected(original,corrected):
    fig = PlotUtilities.figure((3.5,4))
    cmap = plt.cm.afmhot
    ax = plt.subplot(2,1,1)
    imshow_kwargs_original = dict(cmap=cmap,**realistic_min_max(original))
    im = ImageUtil.make_image_plot(original,pct=50,
                                   imshow_kwargs=imshow_kwargs_original)
    ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig)
    PlotUtilities.no_x_label(ax=ax)
    PlotUtilities.xlabel("",ax=ax)
    ax = plt.subplot(2,1,2)
    imshow_kwargs_corrected = dict(cmap=cmap,**realistic_min_max(corrected))
    
    im = ImageUtil.make_image_plot(corrected,pct=50,
                                   imshow_kwargs=imshow_kwargs_corrected)
    ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig,add_space_only=True)  
    return fig
    
def correct_background(images):
    ex_corrected = []
    for ex in images:
        tmp = copy.deepcopy(ex)
        tmp.height = subtract_background(tmp.height)
        ex_corrected.append(tmp)
    return ex_corrected
    
def cache_images(cache_dir,func,**kw):
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
    GenUtilities.ensureDirExists(out_dir)
    images = ImageUtil.cache_images_in_directory(pxp_dir=in_dir,
                                                 cache_dir=cache_dir_raw,
                                                 limit=10)
    corrected_dir = cache_dir_fmt.format("corrected")                                  
    corrected_images = cache_images(corrected_dir,
                                    func = lambda: correct_background(images))
    # subtract the linear backround from each, save to a new cache 
    for i,(orig,corr) in enumerate(zip(images,corrected_images)):
        fig = plot_with_background_corrected(orig,corr)
        img_name = FEC_Util.name_func(i,corr)
        out_name = "{:s}{:s}.png".format(corrected_dir,img_name)
        PlotUtilities.savefig(fig,out_name)

    
    
                                                 
if __name__ == "__main__":
    run()
