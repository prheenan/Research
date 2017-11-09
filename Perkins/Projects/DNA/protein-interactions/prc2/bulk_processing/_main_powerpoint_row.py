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
from GeneralUtil.python import CheckpointUtilities,PlotUtilities,GenUtilities
import ProcessingUtil

def make_row_image(images,out_name,numbers_fmt):
    n_matches = [0 for n in numbers_fmt]
    to_plot = []
    for im in images:
        # see if our name matched the expected format...
        matched = [n in im.Name() for n in numbers_fmt]
        if (matched == [False for _ in matched]):
            continue
        # POST: at least one worked
        for i,m in enumerate(matched):
            n_matches[i] = n_matches[i] + 1 if m else n_matches[i]
        to_plot.append(im)
    # make sure we matched exactly one image
    assert (np.array(n_matches) == 1).all() 
    # plot them all...
    n_cols = len(numbers_fmt)
    size_per_col = 3
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
    to_save = [ ["0x_row.jpeg", [134,144,147,158]],
                ["5mer.jpeg",[50,57,62,110]]]
    for name,numbers in to_save:
        numbers_fmt = ["Image{:04d}".format(n) for n in numbers]
        make_row_image(images,out_dir + name,numbers_fmt)   
    
if __name__ == "__main__":
    run()
