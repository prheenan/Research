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
sys.path.append("../../../../../../../../")
from GeneralUtil.python import CheckpointUtilities

class fig1d:
    def __init__(self,x,y,wlc_x,wlc_y):
        self.x = x
        self.y = y
        self.wlc_x = wlc_x
        self.wlc_y = wlc_y
      
class fig4ab:
    def __init__(self,force,time):
        self.time = time   
        self.force =force
        
class fig4c:
    def __init__(self,energy,energy_error,x):
        self.x = x
        self.energy = energy
        self.energy_error = energy_error 
        
def read_data(input_file,col_func,skip_header=1,delimiter=","):
    cols = np.genfromtxt(base + input_file,skip_header=skip_header,
                         delimiter=delimiter)    
    obj = col_func(cols)
    return obj    
        
def save_output(base,input_file,col_func,**kw):
    output_file = base + input_file + ".pkl"
    obj = CheckpointUtilities.getCheckpoint(output_file,read_data,False,
                                            input_file,col_func,**kw)
    return obj
    
def convert_fig1d(cols):
    n = 8
    all_cols = [col for col in cols.T]
    get_cols = lambda i: all_cols[i*n:(i+1)*n]
    x = get_cols(0)
    y = get_cols(1)
    wlc_y = all_cols[3*n:5*n:2]
    wlc_x = all_cols[3*n+1:5*n:2]
    return fig1d(x=x,y=y,wlc_x=wlc_x,wlc_y=wlc_y)

def run():
    """
    """
    base = "./recreation_figure_data/"
    # # load in 1D ...
    obj = save_output(base=base,input_file="Fig1D.csv",
                      col_func = convert_fig1d)
    # # load in figure 4A,base
    # only load the third and fourth columns (2/3 index)
    func_4 = lambda cols: fig4ab( *([col for col in cols.T][2:4]))
    obj = save_output(base=base,input_file="Fig4AB.csv",
                      col_func = func_4)
    # # load in figure 4C
    obj = save_output(base=base,input_file="Fig4C.csv",
                      col_func = lambda cols: fig4c(*[c for c in cols.T]))

if __name__ == "__main__":
    run()
