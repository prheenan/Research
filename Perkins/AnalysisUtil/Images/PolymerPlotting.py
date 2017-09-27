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

def plot_single_worm_object(tmp_fit):
    plt.subplot(2,1,1)
    plt.plot(*tmp_fit.fit_xy,color='g',marker='o',linewidth=0.5)
    plt.subplot(2,1,2)
    x,y = tmp_fit.x_y_rel
    plt.plot(x,y,color='r',linewidth=2)
    plt.imshow(tmp_fit.image_threshold)
