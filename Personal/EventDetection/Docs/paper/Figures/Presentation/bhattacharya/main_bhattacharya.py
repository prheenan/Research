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

def plot_bhattacharya(sigma,n_samples,bins,loc):
    np.random.seed(42)
    distribution = np.random.normal(loc=0,scale=sigma,size=n_samples)
    uniform = np.random.normal(loc=loc, scale=sigma,size=n_samples)
    min_v,max_v = np.min([min(distribution),min(uniform)]),\
                  np.max([max(distribution),max(uniform)])
    common_style = dict(alpha=0.3,normed=True)
    label_gauss= r"$\varphi$"
    label_uniform=r"$\mathcal{N}$" + \
                 (r"($\mu$={:.1g},$\sigma$={:d})").format(loc,sigma)
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
    loc_arr = [-3,-1,0.5]
    n_bins = 50
    ylim = [0,0.5]
    n_cols = 3+1
    xlim = [-12,12]
    bins = np.linspace(*xlim,endpoint=True,num=n_bins)
    fig = PlotUtilities.figure((16,7))
    plt.subplot(1,n_cols,1)
    bhattacharya = plot_bhattacharya(sigma,n_samples,bins=bins,
                                     loc=-9)
    plt.xlim(xlim)
    PlotUtilities.tickAxisFont()
    PlotUtilities.xlabel("Distribution Value")
    PlotUtilities.ylabel("Probability")
    PlotUtilities.legend(frameon=False)
    plt.ylim(ylim)
    title = (r"BCC of {:.3f}".format(1-bhattacharya))
    PlotUtilities.title(title)
    for i,loc in enumerate(loc_arr):
        plt.subplot(1,n_cols,(i+2))
        bhattacharya = plot_bhattacharya(sigma,n_samples,bins,
                                         loc=loc)
        title = (r"BCC of {:.3f}".format(1-bhattacharya))
        PlotUtilities.title(title)
        plt.xlim(xlim)
        plt.ylim(ylim)
        PlotUtilities.tickAxisFont()
        PlotUtilities.no_y_label()
        PlotUtilities.legend(frameon=False)
        PlotUtilities.xlabel("")
    PlotUtilities.savefig(fig,"bcc.pdf",subplots_adjust=dict(wspace=0.1))
    
if __name__ == "__main__":
    run()
