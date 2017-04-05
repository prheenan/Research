# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os,string

sys.path.append("../../../../../")
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
import cProfile
from Research.Personal.EventDetection._2SplineEventDetector import Detector
def run(base="./"):
    """
    
    """
    name = "examples.pdf"
    data_base = base + "data/"
    file_names = ["fast_unfolding"]
    kw = dict(cache_directory=data_base,force=False)
    file_paths = [data_base + f +".csv" for f in file_names]
    cases = [read_and_cache_file(f,**kw) for f in file_paths]
    ex = cases[0]
    repeats = 10
    pr = cProfile.Profile()
    pr.enable()
    for i in range(repeats):
        _ = Detector.predict(ex,threshold=1e-2)
    pr.disable()
    pr.print_stats(sort='time')

if __name__ == "__main__":
    run()
