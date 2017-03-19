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
    n_sigma = 2
    n_bins = 50
    distribution = np.random.normal(loc=0,scale=sigma,size=n_samples)
    uniform = np.random.uniform(-n_sigma*sigma, n_sigma*sigma, size=n_samples)
    min_v,max_v = np.min([min(distribution),min(uniform)]),\
                  np.max([max(distribution),max(uniform)])
    bins = np.linspace(min_v,max_v,endpoint=True,num=n_bins)
    common_style = dict(alpha=0.3)
    fig = plt.figure(dpi=400)
    n_uniform,e_uniform,_ = plt.hist(uniform,bins=bins,label="uniform",
                                     **common_style)
    n_gauss,e_gauss,_ = plt.hist(distribution,bins=bins,label="gaussian",
                                 hatch="//",**common_style)
    p_uniform = (n_uniform/n_samples)
    p_gauss = (n_gauss/n_samples)
    bhattacharya =  sum(np.sqrt(p_uniform) * np.sqrt(p_gauss))
    title = (r"Uniform over $\pm$" +str(n_sigma) + \
             "$\sigma_{\mathrm{Gaussian}}$"+ \
             " has a BC of {:.3f} with the Gaussian".format(bhattacharya))
    xlim = [min_v,max_v]
    plt.xlim(xlim)
    plt.xlabel("Value")
    plt.ylabel("Count")
    plt.title(title)
    fig.savefig("bc.png")
    
if __name__ == "__main__":
    run()
