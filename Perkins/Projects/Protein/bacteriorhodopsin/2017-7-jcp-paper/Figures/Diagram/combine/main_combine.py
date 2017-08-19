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

import svgutils.transform as sg
from svgutils import compose
import sys

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    f1 = "./diagram.svg"
    f2 = "./Figure3_12_prh.svg"
    size_inches = [7,3]
    size_cm = np.array(size_inches) * 2.54 
    size_x,size_y = ["{:.1f}cm".format(s) for s in size_cm]
    y_text =  sg.TextElement(x=200,y=0.5,text="foo",size=10,
                             weight='bold',font='Arial')
    y_text.rotate(-90,200,0.5)                             
    figure = [compose.Panel(compose.SVG(f1)),
              compose.Panel(compose.SVG(f2)),
              y_text]
    fig = compose.Figure(size_x,size_y,*figure).tile(1,2)
    fig.save('combined.svg')
    
if __name__ == "__main__":
    run()
