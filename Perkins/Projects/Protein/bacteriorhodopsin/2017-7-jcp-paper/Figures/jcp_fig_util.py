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
from GeneralUtil.python.Plot import Annotations
         
def regions_and_colors():
    adhesion_min = 17
    ed_max = 32
    cd_max = 48
    a_max = 65
    regions_colors = [ [[adhesion_min,ed_max],'royalblue'],
                       [[ed_max,cd_max],'orangered'],
                       [[cd_max,a_max],'g']]
    return regions_colors                       
                   
                                  
def add_helical_boxes(ax,alpha=0.3,ymax_box=0.15,font_color=None,
                      offset_bool=False):
    labels_helical_region = ["ED","CB","A"]
    ymin_box = ymax_box-0.1
    regions_colors = regions_and_colors()
    x_arr = [r[0] for r in regions_colors]
    for i,(x,color) in enumerate(regions_colors):
        if (offset_bool):
            x = np.array(x) - np.min(x_arr)
        ax.axvspan(*x,ymin=ymin_box,ymax=ymax_box,color=color,alpha=alpha,
                    linewidth=0)
        ymin,ymax = plt.ylim()
        y_f = (ymin_box+ymax_box)/2 
        y = y_f * (ymax-ymin) + ymin
        x_text = np.mean(x)
        s = labels_helical_region[i]
        font_color_tmp = font_color if font_color is not None else color
        Annotations._annotate(ax=ax,s=s,xy=(x_text,y),
                              horizontalalignment='center',
                              verticalalignment='center',color=font_color_tmp,
                              bbox=dict(alpha=0,pad=0))       