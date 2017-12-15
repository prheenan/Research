# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

class Loops:
    def __init__(self,u,tck,u_starts,u_ends):
        self.u = u
        self.tck = tck
        self.u_starts = u_starts
        self.u_ends = u_ends
    def get_loop_bounds(self):
        return [ [i,f] for i,f in zip(self.u_starts,self.u_ends)]
