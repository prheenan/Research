#!/usr/bin/env python

from libhooke import WX_GOOD
import wxversion
wxversion.select(WX_GOOD)
from wx import PostEvent

import libhookecurve as lhc
from numpy import arange, mean

class superimposeCommands:
    
    def _plug_init(self):
        self.imposed=[]
    
    def do_selimpose(self,args):
        '''
        SELIMPOSE (superimpose.py plugin)
        Hand-selects the curve portion to superimpose
        '''
        #fixme: set to superimpose should be in args
        
        if args=='clear':
            self.imposed=[]
            return
        
        current_set=1
        
        points=self._measure_two_points()
        boundaries=[points[0].index, points[1].index]
        boundaries.sort()
        
        theplot=self.plots[0]
        #append the selected section
        self.imposed.append([])
        self.imposed[-1].append(theplot.vectors[1][0][boundaries[0]:boundaries[1]]) #x
        self.imposed[-1].append(theplot.vectors[1][1][boundaries[0]:boundaries[1]]) #y
        
        #align X first point
        self.imposed[-1][0] = [item-self.imposed[-1][0][0] for item in self.imposed[-1][0]]
        #align Y first point
        self.imposed[-1][1] = [item-self.imposed[-1][1][0] for item in self.imposed[-1][1]]
                
    def do_plotimpose(self,args):
        '''
        PLOTIMPOSE (sumperimpose.py plugin)
        plots superimposed curves
        '''
        imposed_object=lhc.PlotObject()
        imposed_object.vectors=self.imposed
        print 'Plotting',len(imposed_object.vectors),'imposed curves'
        
        imposed_object.normalize_vectors()
        
        imposed_object.units=self.plots[0].units
        imposed_object.title='Imposed curves'
        imposed_object.destination=1
        
        plot_graph=self.list_of_events['plot_graph']       
        PostEvent(self.frame,plot_graph(plots=[imposed_object]))
        
    def do_plotavgimpose(self,args):
        '''
        PLOTAVGIMPOSE (superimpose.py plugin)
        Plots the average of superimposed curves using a running window
        '''
        step=(-5*(10**-10))
        #find extension of each superimposed curve
        min_x=[]
        for curve in self.imposed:
            min_x.append(min(curve[0]))
            
        #find minimum extension
        min_ext_limit=max(min_x)
        
        x_avg=arange(step,min_ext_limit,step)
        y_avg=[]
        for value in x_avg:
            to_avg=[]
            for curve in self.imposed:
                for xvalue, yvalue in zip(curve[0],curve[1]):
                    if xvalue >= (value+step) and xvalue <= (value-step):
                        to_avg.append(yvalue)                                 
            y_avg.append(mean(to_avg))
        
        print 'len x len y'
        print len(x_avg), len(y_avg)
        print y_avg
            
        avg_object=lhc.PlotObject()
        avg_object.vectors=[[x_avg, y_avg]]
        avg_object.normalize_vectors()
        avg_object.units=self.plots[0].units
        avg_object.title="Average curve"
        avg_object.destination=1
         
        plot_graph=self.list_of_events['plot_graph']       
        PostEvent(self.frame,plot_graph(plots=[avg_object]))