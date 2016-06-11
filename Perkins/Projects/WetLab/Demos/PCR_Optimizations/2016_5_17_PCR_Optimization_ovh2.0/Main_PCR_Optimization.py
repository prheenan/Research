# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


sys.path.append("../../../")
from Util import DilutionUtil
from PCR.Optimization.PCR_Opt_Analysis import PCR_Analyze_Objects

class Gradient:
    MACHINE_TOUCHGENE = 0
    MACHINE_BIORAD = 1
    def __init__(self,Temperatures,ConcentrationYields,ConcentrationVolumes,
                 NumRepeats,Description,Machine,Date,NumVialsCombined=1):
        self.Temperatures = np.array(Temperatures)
        self.Concentrations = np.array(ConcentrationYields)
        self.Volumes = np.array(ConcentrationVolumes)
        self.Repeats = NumRepeats
        self.Description = Description
        self.Machine = Machine
        self.Date = Date
        self.NumVialsCombined = NumVialsCombined
    def GetMachineName(self):
        return "TouchGene" if self.Machine==Gradient.MACHINE_TOUCHGENE\
            else "BioRad"
    def GetYieldPer100uLTube(self):
        return self.Volumes * self.Concentrations / self.NumVialsCombined

def run():
    """
    Taken from notebook #2, pp 14
    """
    GradientsObjs = \
        [Gradient(Temperatures=[60.2,62.5,64.2,65.8],
                  ConcentrationYields=[76,62,40,10],
                  ConcentrationVolumes=35,
                  NumRepeats=30,
                  Description="Notebook#2, pp14",
                  Machine=Gradient.MACHINE_TOUCHGENE,
                  Date="???"),
         Gradient(Temperatures=[60.2,62.5,64.2,65.8],
                  ConcentrationYields=[110,100,70,30],
                  ConcentrationVolumes=35,
                  NumRepeats=35,
                  Description="Notebook#2, pp14",
                  Machine=Gradient.MACHINE_TOUCHGENE,
                  Date="???"),
         ## 5/17/2016 data
         # 35R
         Gradient(Temperatures=[60,61.4,62.3,64],
                  ConcentrationYields=[104.1,95.1,96.7,75.7],
                  ConcentrationVolumes=35,
                  NumRepeats=35,
                  Description="",
                  Machine=Gradient.MACHINE_TOUCHGENE,
                  Date="5/17/2016"),
         #40R
         Gradient(Temperatures=[60,61.4,62.3,64],
                  ConcentrationYields=[146.6,149.3,147.4,106.1],
                  ConcentrationVolumes=35,
                  NumRepeats=40,
                  Description="",
                  Machine=Gradient.MACHINE_TOUCHGENE,
                  Date="5/17/2016"),
         # 5/23 data, only one temperature and pooeled two vials
         Gradient(Temperatures=[61.4],
                  ConcentrationYields=[464],
                  ConcentrationVolumes=35,
                  NumRepeats=40,
                  Description="",
                  Machine=Gradient.MACHINE_TOUCHGENE,
                  Date="5/23/2016",
                  NumVialsCombined=2),
         ## 5/24 data, ibid
         Gradient(Temperatures=[61.4],
                  ConcentrationYields=[350],
                  ConcentrationVolumes=35,
                  NumRepeats=40,
                  Description="",
                  Machine=Gradient.MACHINE_TOUCHGENE,
                  Date="5/24/2016",
                  NumVialsCombined=2),
         ## 5/26 data on the biorad
         Gradient(Temperatures=[60,61.3,62.5,64],
                  ConcentrationYields=[155.2,86.7,41.5,50],
                  ConcentrationVolumes=35,
                  NumRepeats=35,
                  Description="",
                  Machine=Gradient.MACHINE_BIORAD,
                  Date="5/26/2016"),
         Gradient(Temperatures=[60,61.3,62.5,64],
                  ConcentrationYields=[172,137,127,62.6],
                  ConcentrationVolumes=35,
                  NumRepeats=40,
                  Description="",
                  Machine=Gradient.MACHINE_BIORAD,
                  Date="5/26/2016"),
         Gradient(Temperatures=[60,61.3],
                  ConcentrationYields=[55,44],
                  ConcentrationVolumes=35,
                  NumRepeats=40,
                  Description="Gradient of ovh-2.0,labelled. T_ann too high?",
                  Machine=Gradient.MACHINE_BIORAD,
                  Date="6/2/2016"),
         Gradient(Temperatures=[58,60.5,62],
                  ConcentrationYields=[95,91,91],
                  ConcentrationVolumes=35,
                  NumRepeats=40,
                  Description="Gradient of ovh-2.0,labelled, with 45s ext",
                  Machine=Gradient.MACHINE_BIORAD,
                  Date="6/3/2016"),
         Gradient(Temperatures=[58,60.5,62],
                  ConcentrationYields=[145,105,121],
                  ConcentrationVolumes=35,
                  NumRepeats=40,
                  Description="Gradient of ovh-2.0,spacer, with 45s ext",
                  Machine=Gradient.MACHINE_BIORAD,
                  Date="6/3/2016"),
         Gradient(Temperatures=[60],
                  ConcentrationYields=[91.5],
                  ConcentrationVolumes=70*4, # 8 tubes into 4, diluted 2-fold
                  NumRepeats=40,
                  Description="Gradient of ovh-2.0,spacer, with 45s ext",
                  Machine=Gradient.MACHINE_BIORAD,
                  Date="6/6/2016")
         
         ]
    PCR_Analyze_Objects(GradientsObjs,"Ovh2.0-Spacer")

if __name__ == "__main__":
    run()
