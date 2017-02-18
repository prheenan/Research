#!/usr/bin/env python

'''
PROCPLOTS
Processed plots plugin for force curves.

Licensed under the GNU GPL version 2
'''
from libhooke import WX_GOOD
import wxversion
wxversion.select(WX_GOOD)

import wx
import libhookecurve as lhc
import numpy as np
import scipy as sp
import scipy.signal
import copy

class procplotsCommands:
    
    def _plug_init(self):
        pass
               
    def do_derivplot(self,args):
        '''
        DERIVPLOT
        (procplots.py plugin)
        Plots the derivate (actually, the discrete differentiation) of the current force curve
        ---------
        Syntax: derivplot
        '''                
        dplot=self.derivplot_curves()               
        plot_graph=self.list_of_events['plot_graph']       
        wx.PostEvent(self.frame,plot_graph(plots=[dplot]))
        
    def derivplot_curves(self):
        '''
        do_derivplot helper function
        
        create derivate plot curves for force curves.
        '''    
        dplot=lhc.PlotObject()
        dplot.vectors=[]       
        
        for vector in self.plots[0].vectors:
            dplot.vectors.append([])
            dplot.vectors[-1].append(vector[0][:-1])
            dplot.vectors[-1].append(np.diff(vector[1]))
        
        dplot.destination=1
        dplot.units=self.plots[0].units
        
        return dplot        
            
    def do_subtplot(self,args):
        '''
        SUBTPLOT
        (procplots.py plugin)
        Plots the difference between ret and ext current curve
        -------
        Syntax: subtplot
        '''
        #FIXME: sub_filter and sub_order must be args
        
        if len(self.plots[0].vectors) != 2:
            print 'This command only works on a curve with two different plots.'
            pass
        
        outplot=self.subtract_curves(sub_order=1)
        
        plot_graph=self.list_of_events['plot_graph']       
        wx.PostEvent(self.frame,plot_graph(plots=[outplot]))
        
    def subtract_curves(self, sub_order):
        '''
        subtracts the extension from the retraction
        ''' 
        xext=self.plots[0].vectors[0][0]
        yext=self.plots[0].vectors[0][1]
        xret=self.plots[0].vectors[1][0]
        yret=self.plots[0].vectors[1][1]
        
        #we want the same number of points
        maxpoints_tot=min(len(xext),len(xret))
        xext=xext[0:maxpoints_tot]
        yext=yext[0:maxpoints_tot]
        xret=xret[0:maxpoints_tot]
        yret=yret[0:maxpoints_tot]
    
        if sub_order:
            ydiff=[yretval-yextval for yretval,yextval in zip(yret,yext)]
        else: #reverse subtraction (not sure it's useful, but...)
            ydiff=[yextval-yretval for yextval,yretval in zip(yext,yret)]    
        
        outplot=copy.deepcopy(self.plots[0])
        outplot.vectors[0][0], outplot.vectors[1][0] = xext,xret #FIXME: if I use xret, it is not correct!
        outplot.vectors[1][1]=ydiff
        outplot.vectors[0][1]=[0 for item in outplot.vectors[1][0]]
        
        return outplot


#-----PLOT MANIPULATORS
    def plotmanip_median(self, plot, current, customvalue=None):
        '''
        does the median of the y values of a plot
        '''
        if customvalue:
            median_filter=customvalue
        else:
            median_filter=self.config['medfilt']
         
        if median_filter==0:
            return plot
         
        if float(median_filter)/2 == int(median_filter)/2:
            median_filter+=1
            
        nplots=len(plot.vectors)
        c=0
        while c<nplots:
            plot.vectors[c][1]=scipy.signal.medfilt(plot.vectors[c][1],median_filter)
            c+=1
        
        return plot
    

    def plotmanip_correct(self, plot, current, customvalue=None):
        '''
        does the correction for the deflection for a force spectroscopy curve.
        Assumes that:
        - the current plot has a deflection() method that returns a vector of values
        - the deflection() vector is as long as the X of extension + the X of retraction
        - plot.vectors[0][0] is the X of extension curve
        - plot.vectors[1][0] is the X of retraction curve
        
        FIXME: both this method and the picoforce driver have to be updated, deflection() must return
        a more senseful data structure!
        '''
        #use only for force spectroscopy experiments!
        if current.curve.experiment != 'smfs':
            return plot
    
        if customvalue != None:
            execute_me=customvalue
        else:
            execute_me=self.config['correct']
        if not execute_me:
            return plot
     
        defl_ext,defl_ret=current.curve.deflection()
        #halflen=len(deflall)/2
    
        plot.vectors[0][0]=[(zpoint-deflpoint) for zpoint,deflpoint in zip(plot.vectors[0][0],defl_ext)]
        plot.vectors[1][0]=[(zpoint-deflpoint) for zpoint,deflpoint in zip(plot.vectors[1][0],defl_ret)]

        return plot


    def plotmanip_centerzero(self, plot, current, customvalue=None):
        '''
        Centers the force curve so the median (the free level) corresponds to 0 N
        Assumes that:
        - plot.vectors[0][1] is the Y of extension curve
        - plot.vectors[1][1] is the Y of retraction curve
        
       
        '''
        #use only for force spectroscopy experiments!
        if current.curve.experiment != 'smfs':
            return plot
    
        if customvalue != None:
            execute_me=customvalue
        else:
            execute_me=self.config['centerzero']
        if not execute_me:
            return plot
     
        
	
	#levelapp=float(np.median(plot.vectors[0][1]))
	levelret=float(np.median(plot.vectors[1][1][-300:-1]))

	level=levelret	

	approach=[i-level for i in plot.vectors[0][1]]
	retract=[i-level for i in plot.vectors[1][1]]
	
	plot.vectors[0][1]=approach	
	plot.vectors[1][1]=retract	
        return plot
    
    '''
    def plotmanip_detriggerize(self, plot, current, customvalue=None):
        #DEPRECATED
        if self.config['detrigger']==0:
            return plot
        
        cutindex=2
        startvalue=plot.vectors[0][0][0]
        
        for index in range(len(plot.vectors[0][0])-1,2,-2):  
           if plot.vectors[0][0][index]>startvalue:
                cutindex=index
           else:
                break

        plot.vectors[0][0]=plot.vectors[0][0][:cutindex]
        plot.vectors[0][1]=plot.vectors[0][1][:cutindex]
        
        return plot
    '''       
           
           
    
#FFT---------------------------
    def fft_plot(self, vector):
        '''
        calculates the fast Fourier transform for the selected vector in the plot
        '''
        fftplot=lhc.PlotObject()
        fftplot.vectors=[[]]
               
        fftlen=len(vector)/2 #need just 1/2 of length
        fftplot.vectors[-1].append(np.arange(1,fftlen).tolist()) 
        
        try:
            fftplot.vectors[-1].append(abs(np.fft(vector)[1:fftlen]).tolist())
        except TypeError: #we take care of newer NumPy (1.0.x)
            fftplot.vectors[-1].append(abs(np.fft.fft(vector)[1:fftlen]).tolist())
        
        
        fftplot.destination=1
        
        
        return fftplot
        
    
    def do_fft(self,args):
        '''
        FFT
        (procplots.py plugin)
        Plots the fast Fourier transform of the selected plot
        ---
        Syntax: fft [top,bottom] [select] [0,1...]
        
        By default, fft performs the Fourier transform on all the 0-th data set on the
        top plot.
        
        [top,bottom]: which plot is the data set to fft (default: top)
        [select]: you pick up two points on the plot and fft only the segment between
        [0,1,...]: which data set on the selected plot you want to fft (default: 0)
        '''
        
        #args parsing
        #whatplot = plot to fft
        #whatset = set to fft in the plot
        select=('select' in args)
        if 'top' in args:
            whatplot=0
        elif 'bottom' in args:
            whatplot=1
        else:
            whatplot=0
        whatset=0
        for arg in args:
            try:
                whatset=int(arg)
            except ValueError:
                pass
        
        if select:
            points=self._measure_N_points(N=2, whatset=whatset)
            boundaries=[points[0].index, points[1].index]
            boundaries.sort()
            y_to_fft=self.plots[whatplot].vectors[whatset][1][boundaries[0]:boundaries[1]] #y
        else:
            y_to_fft=self.plots[whatplot].vectors[whatset][1] #y
        
        fftplot=self.fft_plot(y_to_fft)               
        fftplot.units=['frequency', 'power']
        plot_graph=self.list_of_events['plot_graph']       
        wx.PostEvent(self.frame,plot_graph(plots=[fftplot]))
        
