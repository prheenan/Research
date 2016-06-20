# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../../../")
from Research.Perkins.AnalysisUtil.Gels.ImageJUtil import ReadFileToLaneObj
import GeneralUtil.python.PlotUtilities as pPlotUtil

def run():
    """
    Runs analysis on the pre and post freeze and squeeze for ovh2.0, showing
    that the linear band is mostly linear
    """
    BaseDir = "./Data/"
    LinearPreFile = "Lane3_Linear_Only.xls"
    LinearPostFile = "Lane4_Linear_Post.xls"
    CircularPreFile = "Lane5_Circular_Pre.xls"
    CircularPostFile = "Lane6_Circular_Post.xls"
    LinearPreObj = ReadFileToLaneObj(BaseDir + LinearPreFile)
    LinearPostObj = ReadFileToLaneObj(BaseDir + LinearPostFile)
    CircularPreObj = ReadFileToLaneObj(BaseDir + CircularPreFile)
    CircularPostObj = ReadFileToLaneObj(BaseDir + CircularPostFile)
    # # Get the linear pre and post as a bar plot
    # only index in linear pre is the linear band
    LinearPre = LinearPreObj.Normalized(0)
    # post Linear has both linear and Circular
    print(LinearPostObj)
    LinearPost_Circular = LinearPostObj.Normalized(0)
    LinearPost_Linear = LinearPostObj.Normalized(1)
    # pre circular has concat, circular, linear in that order
    CircularPre_Concat = CircularPreObj.Normalized(0)
    CircularPre_Circular = CircularPreObj.Normalized(2)
    CircularPre_Linear = CircularPreObj.Normalized(2)
    # post circular has concat, circular, linear in that ordr
    CircularPost_Concat = CircularPostObj.Normalized(0)
    CircularPost_Circular = CircularPostObj.Normalized(1)
    CircularPost_Linear = CircularPostObj.Normalized(2)
    # plot them all
    N = 3
    fig = pPlotUtil.figure()
    ax = plt.subplot()
    width = 0.4
    Data = [LinearPre,LinearPost_Circular,LinearPost_Linear,
            CircularPre_Concat,CircularPre_Circular,CircularPre_Linear,
            CircularPost_Concat,CircularPost_Circular,CircularPost_Linear]
    tick_labels = ["PreLinear","PostLinear-Circular","PostLinearLinear",
                   "","","","","",""]
    ax.bar(left=np.arange(len(Data)),
           tick_label=tick_labels,facecolor='yellow',edgecolor='gray')
    pPlotUtil.savefig(fig,"./out.png")
    

if __name__ == "__main__":
    run()
