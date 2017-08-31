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
    f2 = "./Figure3_12_prh_mod.svg"
    size_inches = [7,3]
    size_cm = np.array(size_inches) * 2.54 
    size_x,size_y = ["{:.1f}cm".format(s) for s in size_cm]
    text_label = dict(weight='bold',size=8,font='Arial')
    text_common = dict(weight='bold',size=6,font='Arial')
    text_subplot_label =  dict(weight='bold',size=11,font='Arial')
    red_label = dict(**text_common)
    x_size = 490
    y_size = 300
    y_height = 0.7
    x_height = 1.07
    label_height = 0.52
    y_50 = 0.57
    y_text_arr = [ [0.025,y_height,"Force",text_label],
                   [0.36,y_height,"Force",text_label],
                   [0.69,y_height,"Free Energy",text_label],
                   # scale bars -- y labels 
                   [0.07,y_50,"50 pN",text_common],
                   [0.42,y_50+0.2,"10 pN",text_common],
                   [0.80,y_50+0.22,"0.5 kcal/mol",text_common]
                   ]
    x_text_arr = [ [0.15,x_height,"Extension",text_label],
                   [0.45,x_height,"Time (ms)",text_label],
                   [0.8,x_height,"Extension (nm)",text_label],
                   # plot labels 
                   [0.01,label_height,"D",text_subplot_label],
                   [0.35,label_height,"E",text_subplot_label],
                   [0.68,label_height,"F",text_subplot_label],
                   [0.55,y_50,"10us",red_label],
                   # scale bars -- x labels 
                   [0.08,y_50+0.05,"4nm",text_common],
                   [0.44,y_50+0.23,"50 us",text_common],
                   [0.82,y_50+0.23,"0.1nm",text_common]]
    x_text_elements,y_text_elements = [],[]
    for x,y,s,kw in y_text_arr:
        x *= x_size
        y *= y_size
        y_text = sg.TextElement(x=x,y=y,text=s,**kw)
        y_text_elements.append(y_text)
        y_text.rotate(-90,x,y)      
    for x,y,s,kw in x_text_arr:
        x *= x_size
        y *= y_size
        x_text_elements.append(sg.TextElement(x=x,y=y,text=s,**kw))
    all_elements = x_text_elements + y_text_elements      
    figure = [compose.Panel(compose.SVG(f1)),
              compose.Panel(compose.SVG(f2).scale(0.92)).move(0,y_size*0.47)] + \
              all_elements
    fig = compose.Figure(size_x,size_y,*figure)
    fig.save('combined.svg')
    
if __name__ == "__main__":
    run()
