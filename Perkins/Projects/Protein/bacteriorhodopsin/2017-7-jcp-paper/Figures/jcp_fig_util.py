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
    adhesion_max_nm = 19
    # plot the helical regions...
    ed_end = 29
    cb_end = 45
    regions_colors = [ [[adhesion_max_nm,ed_end],'royalblue'],
                       [[ed_end,cb_end],'orangered'],
                       [[cb_end,67],'g']]
    return regions_colors                       
                   
                                  
def add_helical_boxes(ax,alpha=0.3):
    labels_helical_region = ["ED","CB","A"]
    ymin_box,ymax_box = 0.05,0.15
    regions_colors = regions_and_colors()
    for i,(x,color) in enumerate(regions_colors):
        ax.axvspan(*x,ymin=ymin_box,ymax=ymax_box,color=color,alpha=alpha,
                    linewidth=0)
        ymin,ymax = plt.ylim()
        y_f = (ymin_box+ymax_box)/2 
        y = y_f * (ymax-ymin) + ymin
        x_text = np.mean(x)
        s = labels_helical_region[i]
        Annotations._annotate(ax=ax,s=s,xy=(x_text,y),
                              horizontalalignment='center',
                              verticalalignment='center',color=color,
                              bbox=dict(alpha=0,pad=0))       