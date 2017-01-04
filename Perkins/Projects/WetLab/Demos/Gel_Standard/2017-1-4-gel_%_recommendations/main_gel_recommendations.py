# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../")
from GeneralUtil.python import PlotUtilities

def run():
    """
    plot of how we should run an agarose gel for small fragments

    1.Agilent. Gel Electrophoresis.
    """
    percentage_and_range = [ [0.5,1000,30000],
                             [0.7,800,12000],
                             [1  ,500,10000],
                             [1.2,400,7000],
                             [1.5,200,3000],
                             [2.0,50 ,2000]]
    pct = [p[0] for p in percentage_and_range ]
    bp_low = [p[1] for p in percentage_and_range ]
    bp_high = [p[2] for p in percentage_and_range ]
    # get the coefficients to fit, fitting to a log scale
    coeffs_low = np.polyfit(x=pct,y=np.log(bp_low),deg=1)
    coeffs_high = np.polyfit(x=pct,y=np.log(bp_high),deg=1)
    desired_pct = 2.5
    desired_bp_resolution = 20
    xvals = np.linspace(0,desired_pct)
    # get the fits
    fit_low = np.exp(np.polyval(coeffs_low,xvals))
    fit_high = np.exp(np.polyval(coeffs_high,xvals))
    # get the lower and upper bounds we can resolve at 2.5%
    expected_at_desired_low = np.exp(np.polyval(coeffs_low,desired_pct))
    expected_at_desired_high = np.exp(np.polyval(coeffs_high,desired_pct))
    # make a plot of what we wnat
    fig = PlotUtilities.figure()
    plt.plot(pct,bp_low,'bo-',label="Lower bp data")
    plt.plot(pct,bp_high,'ko-',label="higher bp data")
    plt.plot(xvals,fit_low,'b--',label="Exponential fit")
    plt.plot(xvals,fit_high,'k--')
    plt.plot(desired_pct,expected_at_desired_low,
             color='b',marker="o",markersize=10)
    plt.plot(desired_pct,expected_at_desired_high,
             color='k',marker="o",markersize=10)
    fill_style = dict(alpha=0.3,color='r')
    plt.fill_between(xvals,fit_low,fit_high,**fill_style)
    plt.plot([],[],label="resolvable region",linewidth=20,**fill_style)
    plt.ylim(desired_bp_resolution/2,max(bp_high)*2)
    plt.xlim(0,desired_pct * 1.1)
    plt.yscale('log')
    title = "Expect resolution between {:d}bp and {:d}bp at {:.1f}% Agarose".\
            format(int(expected_at_desired_low),int(expected_at_desired_high),
                   desired_pct)
    xlabel = r"Agarose percentage ($\frac{\mathrm{g}}{\mathrm{mL}}$)"
    PlotUtilities.lazyLabel(xlabel,"Resolvable bp",title,frameon=True)
    PlotUtilities.savefig(fig,"./out.png")

if __name__ == "__main__":
    run()
