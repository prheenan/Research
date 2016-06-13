# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from GeneralUtil.python import CheckpointUtilities as pCheckUtil
from GeneralUtil.python import PlotUtilities as pPlotUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import \
    FEC_Util,FEC_Plot

def FullPlot(Obj):
    FEC_Plot.FEC(Obj,GetX=lambda x: x.Zsnsr,XLabel="Stage Position (nm)",
                 NFilterPoints=40)


def SavePlot(SortByZ,FullName,SaveName):
    Objs = FEC_Util.ReadInData(FullName)
    Tmp = Objs[0]
    if (SortByZ):
        ArgSort = np.argsort(Tmp.Zsnsr)
        Tmp.Zsnsr = Tmp.Zsnsr[ArgSort]
        Tmp.Force = Tmp.Force[ArgSort]
        Offset = 150
    else:
        Offset = 0
    # set up the y limits in pN
    ylim = [Offset-20,Offset+25]
    fig = pPlotUtil.figure()
    plt.subplot(2,1,1)
    FullPlot(Tmp)
    pPlotUtil.lazyLabel("Stage position (nm)","Force (pN)","Flickering for FEC",
                        frameon=True)
    # x limits for 'zoomed in'
    plt.xlim([30,40])
    plt.ylim(ylim)
    plt.subplot(2,1,2)
    FullPlot(Tmp)
    pPlotUtil.lazyLabel("Stage position (nm)","Force (pN)","",frameon=True)
    # x limits for 'zoomed in'
    plt.xlim([36,36.5])
    plt.ylim(ylim)
    pPlotUtil.savefig(fig,SaveName)

    
    
def run():
    Base = "/Users/patrickheenan/Documents/education/boulder_files/" + \
           "rotations_year_1/3_perkins/reports/2016_6_7_flickering_dopc/"
    FullName = Base + "2016_6_7_Experiment_Py.pxp"
    SavePlot(SortByZ=True,FullName=FullName,SaveName = Base + "ZSort.png")
    SavePlot(SortByZ=False,FullName=FullName,SaveName = Base + "Raw.png")
    
if __name__ == "__main__":
    run()
