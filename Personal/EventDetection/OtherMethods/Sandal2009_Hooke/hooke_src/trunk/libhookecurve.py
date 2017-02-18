#!/usr/bin/env python
# -*- coding: utf-8 -*-


class HookeCurve(object):
    
    def __init__(self,path):
        self.path=path
        self.curve=Driver()
        self.notes=''
    
    def identify(self, drivers):
        '''
        identifies a curve and returns the corresponding object
        '''
        for driver in drivers:
	    try:
              tempcurve=driver(self.path)
	    except:
	      print "Error in the playlist of the files."
	      return False
            if tempcurve.is_me():
                #bring on all the driver, with his load of methods etc.
                #so we can access the whole of it.
                self.curve=tempcurve
                del tempcurve
                return True
        
        print 'Not a recognizable curve format.'
        return False
        
        
class Driver:
    '''
    Base class for file format drivers.
    
    To be overridden
    '''
    def __init__(self):
        self.experiment=''
        self.filetype=''
    
    def is_me(self):
        '''
        This method must read the file and return True if the filetype can be managed by the driver, False if not.
        '''
        return False
    
    def close_all(self):
        '''
        This method must close all the open files of the driver, explicitly.
        '''
        return None
    
    def default_plots(self):
        dummy_default=PlotObject()
        dummy_default.vectors.append([[[0]],[[0]]])
        return [dummy_default]
   

class PlotObject:
    
    def __init__(self):
        
        '''
        the plot destination
        0=top
        1=bottom
        '''
        self.destination=0 
        
        '''
        self.vectors is a multidimensional array:
        self.vectors[0]=plot1
        self.vectors[1]=plot2
        self.vectors[2]=plot3
        etc.
        
        2 curves in a x,y plot are:
        [[[x1],[y1]],[[x2],[y2]]]
        for example:
            x1          y1              x2         y2
        [[[1,2,3,4],[10,20,30,40]],[[3,6,9,12],[30,60,90,120]]]
        x1 = self.vectors[0][0]
        y1 = self.vectors[0][1]
        x2 = self.vectors[1][0]
        y2 = self.vectors[1][1]
        '''
        self.vectors=[]

        '''
        self.units is simpler. for each plot with N axes (x,y,z...) only N labels
        can be made, regardless of the number of superimposed plots
        so units for the double plot above is: [unitx, unity]
        
        units are strings
        '''
        self.units=['','']
        
        '''
        xaxes and yaxes directions. 0,0 means the common +X=right, +Y=top directions
        '''
        self.xaxes=0
        self.yaxes=0
        
        self.title='' #title 
        
        '''
        styles: defines what is the style of the current plots. If undefined or None, it is line plot.
        If an element of the list is 'scatter', the corresponding dataset
        is drawn with scattered points and not a continuous line.
        '''
        self.styles=[]
        
        '''
        colors: define what is the colour of the current plots
        '''
        self.colors=[]
        
    def add_set(self,x,y):
        '''
        Adds an x,y data set to the vectors.
        '''
        self.vectors.append([])
        self.vectors[-1].append(x)
        self.vectors[-1].append(y)
        return
    
    def remove_set(self,whichset):
        '''
        Removes a set
        '''
        waste=self.vectors.pop(whichset)
        return
    
    def normalize_vectors(self):
        '''
        Trims the vector lengths as to be equal in a plot.
        '''
        
        for index in range(0,len(self.vectors)):
            vectors_to_plot=self.vectors[index]
            lengths=[len(vector) for vector in vectors_to_plot]
            if min(lengths) != max(lengths):
                for indexplot in range(0,len(vectors_to_plot)):
                    self.vectors[index][indexplot] = self.vectors[index][indexplot][0:min(lengths)]

