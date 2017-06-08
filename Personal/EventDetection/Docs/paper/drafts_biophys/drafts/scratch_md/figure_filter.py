#! /usr/bin/env python3
"""Pandoc filter that replaces labels of format {#?:???}, where ? is a
single lower case character defining the type and ??? is an alphanumeric
label, with numbers. Different types are counted separately.

credit to: blog.hartleygroup.org/2015/11/08/numbering-figures-schemes-and-charts-in-pandoc/

"""

from pandocfilters import toJSONFilter, Str
import re,pickle

REF_PAT = re.compile('(.*)\{#(label_|ref_)?(fig|eq|tbl|sec):(\w*)\}(.*)')

known_labels = {}

def match_pattern(val):
    start, label_or_ref,kind, label, end = REF_PAT.match(val).groups()
    # make the king case-insensitive
    kind = kind.lower()
    return start,label_or_ref,kind,label,end


def determine_references(key, val, fmt, meta):
    if key == 'Str' and REF_PAT.match(val):
        start,label_or_ref,kind,label,end = match_pattern(val)
        if kind in known_labels:
            if label not in known_labels[kind] and (label_or_ref == "label_"):
                known_labels[kind][label] = str(len(known_labels[kind])\
                                                + 1)
        else:
            known_labels[kind] = {}
            known_labels[kind][label] = "1"

def fulfill_references(key,val,fmt,meta):
    if key == 'Str' and REF_PAT.match(val):
        start, label_or_ref,kind, label, end = REF_PAT.match(val).groups()
        return [Str(start)] + [Str(known_labels[kind][label])] + \
               [Str(end)]

if __name__ == '__main__':
    toJSONFilter(determine_references)
    # save the labelling dictionary out as a pkl file
    file_out = "./tmp.pkl"
    with open(file_out,'wb') as f:
        pickle.dump(known_labels,f)
