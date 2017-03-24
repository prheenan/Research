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
sys.path.append("../../../../../../../../")
from GeneralUtil.python import PlotUtilities

def bc(x,y):
    p_x = (x/sum(x))
    p_y = (y/sum(y))
    return sum(np.sqrt(p_x) * np.sqrt(p_y))

def plot_bhattacharya(sigma,n_samples,bins,low,high):
    np.random.seed(42)
    distribution = np.random.normal(loc=0,scale=sigma,size=n_samples)
    uniform = np.random.uniform(low,high, size=n_samples)
    min_v,max_v = np.min([min(distribution),min(uniform)]),\
                  np.max([max(distribution),max(uniform)])
    common_style = dict(alpha=0.3)
    label_gauss= r"$\mathcal{N}$" + \
                 (r"($\mu$={:d},$\sigma$={:d})").format(0,sigma)
    label_uniform=r"$\mathcal{U}$" + (r"(a={:d},b={:d})").format(low,high)
    n_uniform,e_uniform,_ = plt.hist(uniform,bins=bins,label=label_uniform,
                                     **common_style)
    n_gauss,e_gauss,_ = plt.hist(distribution,bins=bins,label=label_gauss,
                                 hatch="////",**common_style)
    return bc(n_uniform,n_gauss)

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    n_samples = int(1e4)
    sigma = 1
    bounds_arr = [[-5,-2],[-3,0],[-2,1]]
    n_bins = 50
    ylim = [0,n_samples/4]
    n_cols = 3+1
    xlim = [-10,10]
    bins = np.linspace(*xlim,endpoint=True,num=n_bins)
    fig = PlotUtilities.figure((16,5))
    plt.subplot(1,n_cols,1)
    bhattacharya = plot_bhattacharya(sigma,n_samples,bins=bins,
                                     low=-9,high=-6)
    plt.xlim(xlim)
    PlotUtilities.tickAxisFont()
    PlotUtilities.xlabel("Value")
    PlotUtilities.ylabel("Count")
    PlotUtilities.legend(frameon=True)
    plt.ylim(ylim)
    title = (r"BC of {:.3f}".format(bhattacharya))
    PlotUtilities.title(title)
    for i,bounds in enumerate(bounds_arr):
        plt.subplot(1,n_cols,(i+2))
        bhattacharya = plot_bhattacharya(sigma,n_samples,bins,
                                         *bounds)
        title = (r"BC of {:.3f}".format(bhattacharya))
        PlotUtilities.title(title)
        plt.xlim(xlim)
        plt.ylim(ylim)
        PlotUtilities.tickAxisFont()
        PlotUtilities.no_y_ticks()
        PlotUtilities.legend(frameon=True)
        PlotUtilities.xlabel("Value")
    PlotUtilities.savefig(fig,"bc.pdf",subplots_adjust=dict(wspace=0.1))
    
if __name__ == "__main__":
    run()
