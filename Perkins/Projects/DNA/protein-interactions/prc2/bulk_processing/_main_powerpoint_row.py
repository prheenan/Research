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
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

from GeneralUtil.python import CheckpointUtilities,PlotUtilities,GenUtilities
import ProcessingUtil

def get_images_to_plot(images,numbers_fmt):
    to_plot = []
    for n in numbers_fmt:
        # see if our name matched the expected format...
        matched = [n in im.Name() for im in images]
        if (matched == [False for _ in matched]):
            continue
        # make sure one image matched
        assert sum(matched) == 1
        # get the matching index
        image_idx = np.where(np.array(matched) == True)[0][0]
        im = images[image_idx]
        to_plot.append(im)
    return to_plot

def make_row_image(to_plot,out_name):
    size_per_col = 3
    n_cols = len(to_plot)
    fig = PlotUtilities.figure((n_cols*size_per_col,3))
    for i,im in enumerate(to_plot):
        is_last = (i == n_cols-1)
        is_first = (i == 0)
        ax = plt.subplot(1,n_cols,(i+1))
        colorbar_kw = dict(add_space_only=not is_last,
                           fontsize=12,fontsize_ticks=12)
        tick_fontsize = 12
        lazy_kwargs = dict(axis_kwargs=dict(fontsize=tick_fontsize))
        ProcessingUtil._image_helper(im,fig,ax,colorbar_kwargs=colorbar_kw,
                                     tick_fontsize=tick_fontsize,
                                     lazy_kwargs=lazy_kwargs)
        if (not is_first):
            PlotUtilities.no_x_label(ax=ax)
            PlotUtilities.no_y_label(ax=ax)
            PlotUtilities.xlabel("",ax=ax)
            PlotUtilities.ylabel("",ax=ax)
        PlotUtilities.tom_ticks(ax=ax,num_major=2)
    subplots_adjust = dict(wspace=-0.02)
    PlotUtilities.savefig(fig,out_name,tight=True,
                          subplots_adjust=subplots_adjust)    

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    images = CheckpointUtilities.lazy_multi_load("./cache/")
    out_dir = "./rows/"
    GenUtilities.ensureDirExists(out_dir)
    preamble = "catalytic-2.5mer"
    to_save = [ ["[prc2]",[79,102,135,97]],
                ["1x",[145,142,132,128]] ,] 
    for name,numbers in to_save:
        numbers_fmt = ["Image{:04d}".format(n) for n in numbers]
        to_plot = get_images_to_plot(images,numbers_fmt)
        base_name = out_dir + preamble + name + ".jpeg"
        make_row_image(to_plot,base_name)
        # save out all the images...
        for i,im in enumerate(to_plot):
            fec_name = FEC_Util.fec_name_func(0,im)
            save_name = "{:s}_{:d}_{:s}".format(base_name,i,fec_name)
            CheckpointUtilities.lazy_save(save_name,im)
            
    
if __name__ == "__main__":
    run()
