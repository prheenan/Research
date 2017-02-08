# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import Bio.SeqUtils
sys.path.append("../../../../../../../")
import GeneralUtil.python.PlotUtilities as pPlotUtil
from Research.Perkins.Projects.WetLab.Util import DepositionUtil

def run_1607F_3520R():
    # get our actual sequence weight
    Primer = "AGAGTGGTCCTA"
    # get the sequence, tack on the overhang
    ProductStartLoc = 1606
    ProductEndLoc = 3520
    with open("mp13_plasmid_plasmid_seq.txt") as f:
        Seq = Primer + "".join([l for l in f.readlines()])
        Seq = Seq[ProductStartLoc:ProductEndLoc]
    # add in the overhang, note I also add an 'A' to represent the abasic site
    Overhang = "TAGGACCACTCT" + "A"
    Seq = Overhang + Seq
    DepositionUtil.run(Seq,output_base="./out")

def run_64nt_hairpin():
    # only add in one side; assume dsDNA
    seq = "gag tca acg tac tga tca cgc tgg atc cta TT"
    # add in the spacers
    seq = "tt" + seq + "tt"
    seq = seq.replace(" ","")
    DepositionUtil.run(seq,output_base="./out_hairpin_")
    
    

def run():
    run_64nt_hairpin()
    run_1607F_3520R()

if __name__ == "__main__":
    run()
