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
from GeneralUtil.python.Plot import Scalebar
import ProcessingUtil

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    images = CheckpointUtilities.lazy_multi_load("./cache_offline/")
    out_dir = "./out_offline/"
    GenUtilities.ensureDirExists(out_dir)
    for i,to_plot in enumerate(images):
        fig = PlotUtilities.figure()
        ax = plt.subplot(1,1,1)
        ProcessingUtil._image_helper(im=to_plot,fig=fig,ax=ax,vmax_nm=2,
                                     vmin_nm=0.1)
        font_kwargs,_ = Scalebar.font_kwargs_modified(dict(color='w'))
        line_kwargs = dict(**Scalebar.def_line_kwargs)
        line_kwargs['color'] = 'w'
        x_unit = PlotUtilities.upright_mu("m")
        Scalebar.x_scale_bar_and_ticks_relative(unit=x_unit,width=1,
                                                offset_x=0.5,
                                                offset_y=0.1,ax=ax,
                                                font_kwargs=font_kwargs,
                                                line_kwargs=line_kwargs)
        PlotUtilities.no_x_anything(ax=ax)
        PlotUtilities.no_y_anything(ax=ax)
        PlotUtilities.savefig(fig,out_dir + str(i) + ".tiff")
            
    
if __name__ == "__main__":
    run()
