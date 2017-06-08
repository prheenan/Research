# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re,pickle

file_cache = "./tmp.pkl"
REF_PAT = re.compile('(.*)\{#(label_|ref_)?((?:S_)?fig|eq|tbl|sec):(\w*)\}(.*)')


def match_pattern(val):
    start, label_or_ref,kind, label, end = REF_PAT.match(val).groups()
    # make the king case-insensitive
    sanit = lambda x : str(x.lower().strip())
    kind = sanit(kind)
    label = sanit(label)
    return start,label_or_ref,kind,label,end

def write_cache(labels):
    with open(file_cache,'wb') as f:
        pickle.dump(labels,f)

def read_cache():
    with open(file_cache,'rb') as f:
        return pickle.load(f)
