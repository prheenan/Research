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


class EnsembleObject:
    def __init__(self, dna_only, dna_plus_protein=None, multimer=None,
                 unknown=None):
        self.dna_only = dna_only
        self.dna_plus_protein = dna_plus_protein
        self.multimer = multimer
        self.unknown = unknown


def convert_to_ensemble(wlc_objs):
    traces = [trace for w in wlc_objs for trace in w.worm_objects]
    dna,dna_plus_loop,multimers,unknown = [],[],[],[]
    for t in traces:
        name = t.file_name.lower()
        if ("protein" in name):
            dna_plus_loop.append(t)
        elif ("dna" in name and "protein" not in name):
            dna.append(t)
        elif ("multiple" in name):
            multimers.append(t)
        elif ("unknown" in name):
            unknown.append(name)
        else:
            assert False , "Didn't recognize {:s}".format(name)
    return EnsembleObject(dna_only=dna,
                          dna_plus_protein=dna_plus_loop,
                          multimer=multimers,
                          unknown=unknown)



