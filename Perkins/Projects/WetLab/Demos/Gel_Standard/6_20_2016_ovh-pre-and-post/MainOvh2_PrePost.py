# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../../../")
from Research.Perkins.AnalysisUtil.Gels.ImageJUtil import \
    ReadFileToOverhangObj,GetLaneTrialsMatchingName
import GeneralUtil.python.PlotUtilities as pPlotUtil

def run():
    """
    Runs analysis on the pre and post freeze and squeeze for ovh2.0, showing
    that the linear band is mostly linear
    """
    BaseDir = "./Data/"
    Files = ["Lane3_Linear_Only.xls","Lane5_Circular_Pre.xls",
             "Lane4_Linear_Post.xls",
             "Lane6_Circular_Post.xls"]
    Labels = ["Linear","Circular DNA, Pre: mixed",
              "Linear Band, Post: mostly linear",
              "Circular Band, Post: mostly Circular",
              "Hot Ligation","Cold Ligation"]
    styles = [dict(color='r',alpha=1),
              dict(color='b',alpha=1.0),
              dict(color='b',alpha=0.6),
              dict(color='b',alpha=0.2)]
    Bins = ["Linear","Circular","Concatemer"]
    Objects = [ReadFileToOverhangObj(BaseDir + f) for f in Files]
    # get the ligaation objects, which have error (multiple trials)
    LigationCold = GetLaneTrialsMatchingName(BaseDir,r""".+Hot_Ligation.+""")
    LigationHot = GetLaneTrialsMatchingName(BaseDir,r""".+Cold_Ligation.+""")
    Objects.extend([LigationHot,LigationCold])
    # construct all the bins
    x_vals = []
    heights = []
    labels = []
    width = 0.4
    # plot them all
    fig = pPlotUtil.figure(figsize=(10,10))
    ax = plt.subplot()
    for i,o in enumerate(Objects):
        intensities = [o.LinearRelative,
                       o.CircularRelative,
                       o.ConcatemerRelative]
        N = len(intensities)
        x = np.arange(0,N) + float(N*i)
        x_vals.extend(x)
        LaneLabel = Labels[i]
        labels.extend([BinLabel
                       for BinLabel in Bins])
        heights.extend(intensities)
        style_idx = i % len(styles)
        try:
            errs = o.GetErrors()
            error_kw = dict(ecolor='k',linewidth=3)
        except AttributeError:
            errs = 0
            error_kw = dict()
        plt.bar(left=x,height=intensities,yerr=errs,error_kw=error_kw,
                edgecolor='gray',label=LaneLabel,
                **(styles[style_idx]))
    span_style = [dict(color='k',alpha=0.3),
                  dict(color='w',alpha=0.0)]
    label_x = np.array( x_vals) + width
    plt.xticks(label_x,labels,rotation=60)
    plt.ylim([0,1.15])
    pPlotUtil.lazyLabel("Population","Population Fraction by Intensity",
                "Purifying circular DNA reveals a dominantly linear population",
                        frameon=True)
    for i,_ in enumerate(Objects):
        start = N*i
        end = N*(i+1)
        color = 'k' if (i % 2) == 0 else "w"
        plt.axvspan(start,end,color=color,alpha=0.2)
    pPlotUtil.savefig(fig,"./out.png")
    

if __name__ == "__main__":
    run()
