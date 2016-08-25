# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
import GeneralUtil.python.PlotUtilities as pPlotUtil

def run():
    """
    """
    DataFile = "./data.txt"
    data = np.loadtxt(DataFile,delimiter=",",dtype=str)
    # load in the character names, removing junk characters
    SkipCharactersForNames = 4
    NameCol = 0
    Names = [d[SkipCharactersForNames:] for d in data[:,NameCol]]
    # load in all scores
    ScoreStartCol = 1
    ScoreEndCol = 5
    # need to convert from table (giving percentage) to just votes
    ConvertStringToScore = lambda x: float(x.split(" ")[1])
    ScoresByNameAndVoteType =[
        [ConvertStringToScore(score)for score in ScoresForThisName]
        for ScoresForThisName in data[:,ScoreStartCol:ScoreEndCol]]
    # add in the weights and sum
    Weights = np.array([5,3,1,-2])
    ScorePerNameToTotal = lambda x: sum(np.array(x) * Weights)
    ScoresByName = [ScorePerNameToTotal(d) for d in ScoresByNameAndVoteType]
    Max = np.max(ScoresByName)
    Min = np.min(ScoresByName)
    # get the y tick lengths
    NumNames= len(Names)
    bins = np.arange(NumNames)
    fig = pPlotUtil.figure(figsize=(8,12))
    ax =plt.subplot(1,1,1)
    IsWinner = lambda x: x == Max
    colors = ['r' if IsWinner(s) else 'k' for s in ScoresByName]
    linewidth = [3 if IsWinner(s) else 0 for s in ScoresByName]
    edgecolor = ['g' if IsWinner(s) else 'k' for s in ScoresByName]
    ax.barh(bins,ScoresByName,align="center",alpha=0.6,color=colors,
            linewidth=linewidth)
    ax.set_yticklabels(Names)
    ax.set_yticks(bins)
    fudge = (Max-Min)/10
    plt.xlim([Min-fudge,Max+fudge])
    plt.axvline(Max,linewidth=3.0,linestyle="--",color="g",
                label="Winning score: {:d}".format(int(Max)))
    pPlotUtil.lazyLabel("Score","Proposed Name Pair","Cypher voting results!")
    pPlotUtil.savefig(fig,"./out.png")
if __name__ == "__main__":
    run()
