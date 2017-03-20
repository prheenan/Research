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
    n_sigma = 5
    distribution = np.random.normal(loc=0,scale=sigma,size=n_samples)
    min_v,max_v = np.min(distribution),np.max(distribution)
    uniform = np.random.uniform(-n_sigma*sigma, n_sigma*sigma, size=n_samples)
    bins = np.linspace(min(uniform),max(uniform),endpoint=True,num=50)
    common_style = dict(alpha=0.3)
    fig = plt.figure(dpi=400)
    n_uniform,e_uniform,_ = plt.hist(uniform,bins=bins,label="uniform",
                                     **common_style)
    n_gauss,e_gauss,_ = plt.hist(distribution,bins=bins,label="gaussian",
                                 hatch="//",**common_style)
    p_uniform = (n_uniform/n_samples)
    p_gauss = (n_gauss/n_samples)
    bhattacharya =  sum(np.sqrt(p_uniform) * np.sqrt(p_gauss))
    title = (r"Uniform over $\pm$5$\sigma_{\mathrm{Gaussian}}$"+ \
             " has a BC of {:.3f} with the Gaussian".format(bhattacharya))
    plt.xlabel("Value")
    plt.ylabel("Count")
    plt.title(title)
    fig.savefig("bc.png")
    
if __name__ == "__main__":
    run()
