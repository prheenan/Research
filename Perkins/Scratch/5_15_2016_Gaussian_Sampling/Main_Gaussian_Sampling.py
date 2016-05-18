# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.stats.distributions import norm



def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    nSamples = [10,50,150,300,400,500,1000]
    distArr = []
    np.random.seed(42)
    for n in nSamples:
        dist = np.random.normal(loc=0.0, scale=1.0, size=n)
        distArr.append(dist)
    # make a figure showing the different distributions for
    # different number of samples
    nBins = 30
    bins = np.linspace(-3,3,nBins)
    nGraphs = len(nSamples)
    fig = plt.figure(figsize=(8,24))
    plt.subplot(nGraphs,1,1)
    plt.title("Recommend at least 300 points from N(0,1)",
              fontsize=25)
    for i,n in enumerate(nSamples):
        plt.subplot(nGraphs,1,i+1)
        plt.hist(distArr[i],bins=bins,label="N={:d}".format(n),normed=True)
        labelPDF = "N(0,1) Dist." if (i == nGraphs-1) else ""
        plt.plot(bins,norm.pdf(bins),label=labelPDF,
                 linewidth=3,linestyle='--',color='r')
        plt.legend()
    plt.ylabel("Weighted Proportion",fontsize=20)
    plt.xlabel("Value of Standard Normal",fontsize=20)
    plt.tight_layout()
    fig.savefig("CompareDist.png")
        

if __name__ == "__main__":
    run()
