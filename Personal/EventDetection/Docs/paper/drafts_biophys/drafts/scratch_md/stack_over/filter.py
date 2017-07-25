
from pandocfilters import toJSONFilter, Str,Emph,Strong

def boldify(key, val, fmt, meta):
    if key == 'Str' and "foo" in val:
        # this is the part I can't don't know how to do
        return [Strong([Str("")])]
        
if __name__ == '__main__':
    toJSONFilter(boldify)
