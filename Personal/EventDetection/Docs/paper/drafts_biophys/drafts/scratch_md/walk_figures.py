#! /usr/bin/env python3
"""Pandoc filter that replaces labels of format {#?:???}, where ? is a
single lower case character defining the type and ??? is an alphanumeric
label, with numbers. Different types are counted separately.

credit to: blog.hartleygroup.org/2015/11/08/numbering-figures-schemes-and-charts-in-pandoc/

"""

from pandocfilters import toJSONFilter, Str
import re,pickle

import figure_util
prefix = [""]
known_labels = {}

def get_prefix():
    return prefix[0]

def set_prefix(s):
    prefix[0] = s


def determine_references(key, val, fmt, meta):
    if key == 'Str' and figure_util.REF_PAT.match(val):
        start,label_or_ref,kind,supp,label,end = figure_util.match_pattern(val)
        if (supp is not None):
            # then we mark all labelled prefixes after as being "S",
            # and remove this meta information
            set_prefix("S")
        # POST: a 'proper' label
        if kind in known_labels:
            if label not in known_labels[kind] and (label_or_ref == "label_"):
                n_labels = len(known_labels[kind])
                known_labels[kind][label] =  get_prefix() + str(n_labels+ 1)
        else:
            known_labels[kind] = {}
            known_labels[kind][label] =  get_prefix() + "1"

if __name__ == '__main__':
    toJSONFilter(determine_references)
    # save the labelling dictionary out as a pkl file
    figure_util.write_cache(known_labels)
