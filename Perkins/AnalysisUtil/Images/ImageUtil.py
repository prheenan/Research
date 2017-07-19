# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


def PlotImage(Image,height_to_plot=None,fix_extent=False,aspect='auto',
              cmap=plt.cm.Greys,range_plot=None,ax=None,**kwargs):
    """
    Plots SurfaceImage as a greyscale map
    
    Args:
         Image:  output of ReadImageAsObject
         label: what to label (title) this with
         **kwargs: passed to imshow
    """
    if (range_plot is None):
        range_plot = Image.range_microns()
    if (height_to_plot is None):
        height_to_plot =Image.height_nm_rel()
    if (fix_extent):
        extent=[0,0,range_plot,range_plot]
    else:
        extent = None
    if (ax is None):
        ax = plt.gca()
    im = ax.imshow(height_to_plot,extent=extent,
                   cmap=cmap,aspect=aspect,**kwargs)
    ax.invert_yaxis()
    # remove the ticks
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right='off', left='off')
    return im

def PlotImageDistribution(Image,pct=95,bins=100,PlotLines=True,AddSigmas=True,
                          **kwargs):
    """
    Plots the distribution of image heights, relative to the surface

    Args:
        Image: output of ReadImageAsObject
        label: see PlotImage
        pct: where we want to draw the line on the distribution, after the
        cummulative probability is over this number [0,100]
        
        bins: number of bins for the histogram
        kwargs: passed to hist
    Returns:
        tuple of <n,bins,patches>, se matplotlib.pyplot.hist
    """
    height_nm_relative = Image.height_nm_rel()
    n,bins,patches = plt.hist(height_nm_relative.ravel(),bins=bins,linewidth=0,
                              edgecolor="none",alpha=0.3,
                              normed=False,**kwargs)
    bin_width = np.median(np.diff(bins))
    height_nm_rel_encompassing_pct = np.percentile(height_nm_relative,pct)
    if (PlotLines):
        min_height = np.min(height_nm_relative.ravel())
I        max_height = np.max(height_nm_relative.ravel())
        limits = np.linspace(start=min_height,stop=max_height,num=bins.size)
        # fit symmetrically to between pct_min% and (100-pct_min)%
        pct_min = 5
        pct_max = 80
        q_min,q_max = np.percentile(height_nm_relative,[pct_min,pct_max])
        fit_idx = np.where( (height_nm_relative > q_min) &
                            (height_nm_relative < q_max))
        height_fit = height_nm_relative[fit_idx]
        mu,std = norm.fit( height_fit)
        n_points = height_fit.size
        denormalize = n_points * bin_width
        pdf = norm.pdf(limits,loc=mu,scale=std) * denormalize
        plt.plot(limits,pdf,label=("Gaussian, stdev={:.2f}".format(std)))
    plt.yscale('log')
    n_limits = n[np.where(n>=1)]
    plt.ylim([min(n_limits)/2,max(n_limits)*2])
    return n,bins,patches

