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
    Labels = ["Linear","Pre-Purification is mixed",
              r"Linear Band is linear",
              r"Circular Band is circular",
              "Ligation, Hot Gel","Ligation, Cold Gel"]
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
    # create a matrix to hold all the measurements
    NumTypes = 3 # linear, circular, concatemer
    NumObjects = len(Objects)
    Matrix = np.zeros((NumObjects,NumTypes))
    for i,o in enumerate(Objects):
        Matrix[i,:] = [o.LinearRelative,
                       o.CircularRelative,
                       o.ConcatemerRelative]
    bins_x = np.arange(NumObjects)
    bins_y = np.arange(NumTypes)
    fig = pPlotUtil.figure(figsize=(14,10))
    ax = plt.subplot(1,1,1)
    heatmap = plt.pcolor(Matrix,cmap=plt.cm.afmhot)
    # put the major ticks at the middle of each cell
    # see:
    # stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor
    ax.set_xticks(np.arange(Matrix.shape[1])+0.5, minor=False)
    ax.set_yticks(np.arange(Matrix.shape[0])+0.5, minor=False)
    # want a more natural, table-like display
    ax.invert_yaxis()
    #ax.xaxis.tick_top()
    # fix the labels
    row_labels = Bins
    column_labels = Labels
    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    cbar = plt.colorbar()
    cbar.set_label('Fraction of population in lane',
                   labelpad=20,rotation=270,size=22)
    cbar.ax.tick_params(labelsize=18)
    title = r"Some DNA remains linear after purification or ligation, " +\
            "suggesting non-annealing DNA"
    pPlotUtil.lazyLabel("DNA Population","Gel Contents",
                        title,frameon=True,titley=1.02)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left='off',      # ticks along the bottom edge are off
        right='off',         # ticks along the top edge are off
        labelleft='on') # labels along the bottom edge are off
    pPlotUtil.savefig(fig,"./out.png")
    

if __name__ == "__main__":
    run()
