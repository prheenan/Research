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


    --According to Figure 2 (pp5) of :

    Mitas, M. et al. 
    Hairpin properties of single-stranded DNA containing a GC-rich triplet 
    repeat: (CTG)15. Nucleic Acids Res 23, 1050-1059 (1995).

    --an 8% polyacrilimide gel works well to resolve ssDNA and dsDNA on the 
    order of 15bp. (we are dealing with more like 30). According to 

    1.Agilent. Gel Electrophoresis

    --an 8% gel is 50-400bp for resolving. The script you are reading
    shows that a 2.5% agarose gel should behave similarly. Note that by ibid,
    the bromophenol blue goes with 45bp (ie: probably want to stop when it is
    <= 1/2 for our ~30bp) 

    See also: 

    Ausubel, F. M. Current protocols in molecular biology, 1987
    """
    percentage_and_range = [ [0.5,1000,30000],
                             [0.7,800,12000],
                             [1  ,500,10000],
                             [1.2,400,7000],
                             [1.5,200,3000],
                             [2.0,50 ,2000],
                             [3.0,10,1000]]
    pct = [p[0] for p in percentage_and_range ]
    bp_low = [p[1] for p in percentage_and_range ]
    bp_high = [p[2] for p in percentage_and_range ]
    # get the coefficients to fit, fitting to a log scale
    coeffs_low = np.polyfit(x=pct,y=np.log(bp_low),deg=1)
    coeffs_high = np.polyfit(x=pct,y=np.log(bp_high),deg=1)
    desired_pct = 2.5
    desired_bp_resolution = 20
    max_pct = max(pct)
    max_graph = max(max_pct,desired_pct)
    xvals = np.linspace(0,max_graph)
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
    plt.xlim(0,max_graph * 1.1)
    plt.yscale('log')
    title = "Expect resolution between {:d}bp and {:d}bp at {:.1f}% Agarose".\
            format(int(expected_at_desired_low),int(expected_at_desired_high),
                   desired_pct)
    xlabel = r"Agarose percentage ($\frac{\mathrm{g}}{\mathrm{mL}}$)"
    PlotUtilities.lazyLabel(xlabel,"Resolvable bp",title,frameon=True)
    PlotUtilities.savefig(fig,"./out_agarose.png")
    """
    # seems a bit sketchy, but here is some data on 
    http://openwetware.org/wiki/Bromophenol_blue
    1000bp at 0.7% agarose
    500bp at 1%
    150bp at 2%
    50bp at 3% 
    """
    pct_bp = [ [0.7,1000],
               [1  ,500 ],
               [2  ,150 ],
               [3  ,50  ]]
    pct = [p[0] for p in pct_bp]
    bp = [p[1] for p in pct_bp]
    lin_coeffs_dye = np.polyfit(pct,np.log(bp),deg=1)
    x_fit = np.linspace(min(pct),max(pct))
    fit_dye = np.exp(np.polyval(lin_coeffs_dye,x_fit))
    expected_running = np.exp(np.polyval(lin_coeffs_dye,desired_pct))
    fig = PlotUtilities.figure()
    plt.plot(pct,bp,'ro-')
    plt.plot(x_fit,fit_dye,'b--')
    plt.plot(desired_pct,expected_running,'ko',markersize=10)
    title = "Effective size (in bp) of Bromophenol blue ({:d} at {:.1f}%)".\
            format(int(expected_running),desired_pct)
    PlotUtilities.lazyLabel("percentage","Effective dye size (bp)",title)
    plt.yscale('log')
    plt.ylim(min(bp)/2,max(bp)*2)
    plt.xlim(0,max_graph*1.1)
    PlotUtilities.savefig(fig,"./out_dye.png")

if __name__ == "__main__":
    run()
