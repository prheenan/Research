#! /usr/bin/env python3
"""Pandoc filter that replaces labels of format {#?:???}, where ? is a
single lower case character defining the type and ??? is an alphanumeric
label, with numbers. Different types are counted separately.

credit to: blog.hartleygroup.org/2015/11/08/numbering-figures-schemes-and-charts-in-pandoc/

"""

from pandocfilters import toJSONFilter, Str
import re

REF_PAT = re.compile('(.*)\{#(fig|eq|tbl|sec):(\w*)\}(.*)')

known_labels = {}

def figref(key, val, fmt, meta):
    if key == 'Str' and REF_PAT.match(val):
        start, kind, label, end = REF_PAT.match(val).groups()
        if kind in known_labels:
            if label not in known_labels[kind]:
                known_labels[kind][label] = str(len(known_labels[kind])\
                                                + 1)
        else:
            known_labels[kind] = {}
            known_labels[kind][label] = "1"
        return [Str(start)] + [Str(known_labels[kind][label])] + \
               [Str(end)]

if __name__ == '__main__':
    toJSONFilter(figref)
