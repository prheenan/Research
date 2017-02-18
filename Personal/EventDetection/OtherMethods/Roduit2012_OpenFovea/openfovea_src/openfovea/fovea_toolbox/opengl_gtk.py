#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
    This module is based on teapot2 example from the python-gtkglext examples.
"""

import pygtk
pygtk.require('2.0')
from gtk.gtkgl.apputils import *

from opengl import Topography, Tomography

class TopoWindow(gtk.Window):
    """
        Displays the 3D array in a gtk window.
        
        Parameters :
            array : 2D numpy.array
                Represents the height
            color_array : 2D numpy.array
                Represents value for the colorscale to plot in the surface.
    """
    def __init__(self, array, color_array=None, clut=None):
        self.key = None
        gtk.Window.__init__(self)
        # Set self attributes.
        self.set_title('3D view')
        if sys.platform != 'win32':
            self.set_resize_mode(gtk.RESIZE_IMMEDIATE)
        self.set_reallocate_redraws(True)
        self.connect('destroy', lambda quit: gtk.main_quit())
        # Catch keyboard events
        self.connect("key-press-event",self.on_window_key_press_event)
        self.connect("key-release-event",self.on_window_key_release_event)
        # Create the table.
        self.table = gtk.Table(3, 3)
        self.table.set_border_width(5)
        self.table.set_col_spacings(5)
        self.table.set_row_spacings(5)
        self.table.show()
        self.add(self.table)
        # The scene and the
        # GLArea widget to
        # display it.
        if array.ndim == 2:
            self.topo = Topography(array, color_array, clut)
        elif array.ndim ==3:
            self.topo = Tomography(array, clut)
        self.glarea = GLArea(self.topo)
        self.glarea.connect("scroll_event", self.button_scroll)
        self.glarea.set_size_request(300,300)
        self.glarea.show()
        self.table.attach(self.glarea, 0, 1, 0, 1)

    def on_window_key_press_event(self, window, event):
        if self.key is None and (64 < event.keyval < 123):
            self.key = gtk.gdk.keyval_name(event.keyval)
    def on_window_key_release_event(self, window, event):
        self.key = None
    def button_scroll(self, window, event):
        if event.direction == gtk.gdk.SCROLL_UP:
            if self.key == None:
                self.topo.zoom = self.topo.zoom + self.topo.zoom/10.
            elif self.key in ['x', 'X', 'y', 'Y', 'z', 'Z']:
                self.topo.slice_it(self.key, -1)
        elif event.direction == gtk.gdk.SCROLL_DOWN:
            if self.key == None:
                self.topo.zoom = self.topo.zoom - self.topo.zoom/10.
            elif self.key in ['x', 'X', 'y', 'Y', 'z', 'Z']:
                self.topo.slice_it(self.key, +1)
        self.glarea.window.invalidate_rect(self.glarea.allocation, False)
    def vchanged(self, vadj):
        #self.topo.rotx = vadj.value
        #self.glarea.window.invalidate_rect(self.glarea.allocation, False)
        self.glarea.window.clear()
        self.topo.slice_it([vadj.value, 0, 0], [-1, -1, -1])
        #help(self.glarea.window.invalidate_rect)
        self.glarea.window.invalidate_rect(self.glarea.allocation, False)
    def hchanged(self, hadj):
        self.topo.roty = hadj.value
        self.glarea.window.invalidate_rect(self.glarea.allocation, False)
    def zoom_changed(self, adj):
        self.topo.zoom = adj.value
        self.glarea.window.invalidate_rect(self.glarea.allocation, False)
    
    def run(self):
        self.show()
        gtk.main()
def main():
    import numpy
    
    size = 32
    # Just create a beatifull gride...
    def func3(x,y):
        return (1- x/2 + x**5 + y**3)*numpy.exp(-x**2-y**2)
    # make these smaller to increase the resolution
    dx, dy = 6./size, 6./size#0.05, 0.05

    x = numpy.arange(-3.0, 3.0, dx)
    y = numpy.arange(-3.0, 3.0, dy)
    X,Y = numpy.meshgrid(x, y)
    
    array = func3(X, Y)
    
    #color_array = numpy.ma.array(array, mask = array<-0.8)
    #app = TopoWindow(array, color_array)
    #app.run()
    array3 = [array*i for i in range(10)]
    array3 = numpy.array(array3).transpose(2,1,0)
    array3 = numpy.ma.array(array3, mask = array3<-0.5)
    clut = (0, 1)
    app = TopoWindow(array3, clut=clut)
    app.run()
if __name__ == '__main__':
    main()
