# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os

sys.path.append("../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from Research.Personal.EventDetection.Util import InputOutput
from GeneralUtil.python import GenUtilities,CheckpointUtilities

            
def run():
    """
    utility process which reads in asylum-style pxp files and converts them and
    their events into csv files. 
    """
    network = FEC_Util.default_data_root()
    base_directory= network + "4Patrick/CuratedData/"
    output_base_directory = base_directory + "Masters_CSCI/"
    positive_directory = output_base_directory + "Positive/"
    relative_650 = "650nm-4x-bio/pxp/"
    relative_protein = "4nug2_alpha3D_MarcAndre/pxp/"
    relative_devin = "4nug2-devin/pxp/"
    relative_input_dir = [relative_protein,
                          relative_devin,
                          relative_650 + "1000-nanometers-per-second/",
                          relative_650 + "100-nanometers-per-second/",
                          relative_650 + "500-nanometers-per-second/"]
    absolute_input_dir = [positive_directory + d 
                          for d in [relative_input_dir[0]]]
    absolute_output_dir = [d.replace("pxp","csv") for d in absolute_input_dir]
    for i,(d,d_out) in enumerate(zip(absolute_input_dir,absolute_output_dir)):
        InputOutput.output_waves_in_directory_to_csv_files(d,d_out)
        break

if __name__ == "__main__":
    run()
