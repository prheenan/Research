#! /usr/bin/env python3
"""Pandoc filter that replaces labels of format {#?:???}, where ? is a
single lower case character defining the type and ??? is an alphanumeric
label, with numbers. Different types are counted separately.

credit to: blog.hartleygroup.org/2015/11/08/numbering-figures-schemes-and-charts-in-pandoc/

"""


from pandocfilters import toJSONFilter, Str,Strong
import re,pickle

import figure_util

def label_text(kind,label,known_labels):
    default_fmt = "{:s} {:s}."
    num = known_labels[kind][label]
    kind_to_label = dict(fig=default_fmt.format("Figure",num),
                         eq=default_fmt.format("Equation",num),
                         tbl=default_fmt.format("Table",num),
                         sec="{:s} {:s}:".format("Section",num))
    return [Strong([Str(kind_to_label[kind])])]
                         

def fulfill_references(key,val,fmt,meta,known_labels):
    if (key == 'Str') and figure_util.REF_PAT.match(val):
        start, label_or_ref,kind, label, end = figure_util.match_pattern(val)
        if (label_or_ref is not None) and ("label" in label_or_ref):
            # format labels as we want...
            content = label_text(kind,label,known_labels)
        elif ( (kind  not in known_labels) or 
             (label not in known_labels[kind])):
            # couldnt find the kind or the label for the reference
            error = "XXX_Unknown-kind/label ({:s}:{:s})XXX".format(kind,label)
            content = [Str(error)]
        else:
            # just a normal reference
            content = [Str(known_labels[kind][label])]
        return [Str(start)] + content + [Str(end)]

if __name__ == '__main__':
    known_labels = figure_util.read_cache()
    with open("tmp.txt",'w') as f:
        f.write(str(known_labels))
    m_func = lambda *args,**kwargs: \
             fulfill_references(*args,
                                known_labels=known_labels,
                                **kwargs)
    toJSONFilter(m_func)
