#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
generalvclamp.py

Plugin regarding general velocity clamp measurements
'''

from libhooke import WX_GOOD, ClickedPoint
import wxversion
wxversion.select(WX_GOOD)
from wx import PostEvent
import numpy as np
import scipy as sp
import copy
import os.path
import time

import warnings
warnings.simplefilter('ignore',np.RankWarning)


class generalvclampCommands:
    
    def _plug_init(self):
        self.basecurrent=None
        self.basepoints=None
        self.autofile=''
    
    def do_distance(self,args):
        '''
        DISTANCE
        (generalvclamp.py)
        Measure the distance (in nm) between two points.
        For a standard experiment this is the delta X distance.
        For a force clamp experiment this is the delta Y distance (actually becomes
        an alias of zpiezo)
        -----------------
        Syntax: distance
        '''
        if self.current.curve.experiment == 'clamp':
            print 'You wanted to use zpiezo perhaps?'
            return
        else:
            dx,unitx,dy,unity=self._delta(set=1)
            print "%.1f nm" %(dx*(10**9))
            to_dump='distance '+self.current.path+' '+str(dx*(10**9))+' nm'
            self.outlet.push(to_dump)


    def do_force(self,args):
        '''
        FORCE
        (generalvclamp.py)
        Measure the force difference (in pN) between two points
        ---------------
        Syntax: force
        '''    
        if self.current.curve.experiment == 'clamp':
            print 'This command makes no sense for a force clamp experiment.'
            return
        dx,unitx,dy,unity=self._delta(set=1)
        print " %.1f pN"%(dy*(10**12))
        to_dump='force '+self.current.path+' '+str(dy*(10**12))+' pN'
        self.outlet.push(to_dump)
        
        
    def do_forcebase(self,args):
        '''
        FORCEBASE
        (generalvclamp.py)
        Measures the difference in force (in pN) between a point and a baseline
        took as the average between two points.
        
        The baseline is fixed once for a given curve and different force measurements,
        unless the user wants it to be recalculated
        ------------
        Syntax: forcebase [rebase]
                rebase: Forces forcebase to ask again the baseline
                max: Instead of asking for a point to measure, asks for two points and use
                     the maximum peak in between
        '''
        rebase=False #if true=we select rebase
        maxpoint=False #if true=we measure the maximum peak
        
        plot=self._get_displayed_plot()
        whatset=1 #fixme: for all sets
        if 'rebase' in args or (self.basecurrent != self.current.path):
            rebase=True
        if 'max' in args:
            maxpoint=True
               
        if rebase:
            print 'Select baseline'
            self.basepoints=self._measure_N_points(N=2, whatset=whatset)
            self.basecurrent=self.current.path
        
        if maxpoint:
            print 'Select two points'
            points=self._measure_N_points(N=2, whatset=whatset)
            boundpoints=[points[0].index, points[1].index]
            boundpoints.sort()
            try:
                y=min(plot.vectors[whatset][1][boundpoints[0]:boundpoints[1]])
            except ValueError:
                print 'Chosen interval not valid. Try picking it again. Did you pick the same point as begin and end of interval?'
        else:
            print 'Select point to measure'
            points=self._measure_N_points(N=1, whatset=whatset)
            #whatplot=points[0].dest
            y=points[0].graph_coords[1]
        
        #fixme: code duplication
        boundaries=[self.basepoints[0].index, self.basepoints[1].index]
        boundaries.sort()
        to_average=plot.vectors[whatset][1][boundaries[0]:boundaries[1]] #y points to average
        
        avg=np.mean(to_average)
        forcebase=abs(y-avg)
        print "%.1f pN"%(forcebase*(10**12))
        to_dump='forcebase '+self.current.path+' '+str(forcebase*(10**12))+' pN'
        self.outlet.push(to_dump)


    def do_aveforce(self,args):
        '''
        AVEFORCE
        (generalvclamp.py)
        Measures the average force for two region and give in output their difference
	with an estimation of the error (using the standard deviation of the average).
        The error is the sum of the standard deviation of the average of the two region.

        ------------
        Syntax: aveforce [set]
                set: Integer number that represent the graph on which we want to operate.
   
        '''

        
        plot=self._get_displayed_plot()
	
	
        whatset=0 #fixme: for all sets
	if len(args)== 0:
	  whatset=0
	if len(args)>0:
	  try:
	    whatset=int(args)
	  except:
	    print "Non integer number in argument."
	    return 0
       
	control=1
	print 'Select baseline'
	while(control==1):
          self.basepoints=self._measure_N_points(N=2, whatset=whatset)
          self.basecurrent=self.current.path
	  if(abs(self.basepoints[0].index-self.basepoints[1].index)>2):
	    control=0
	  else:
	    print "The region is too short, try again."

        control=1
        #we find the two index of the baseline
        boundaries=[self.basepoints[0].index, self.basepoints[1].index]
        boundaries.sort()
        to_base_average=plot.vectors[whatset][1][boundaries[0]:boundaries[1]] #y points to average

	#calculating average and standard deviation of the average for the baseline
        avg_base=np.mean(to_base_average)
	std_avg_base=np.std(to_base_average)/np.sqrt(len(to_base_average))
	

	print 'Select two points'
	while(control==1):
	  points=self._measure_N_points(N=2, whatset=whatset)
	  if(abs(points[0].index-points[1].index)>2):
	    control=0
	  else:
	    print "The region is too short, try again."

        boundpoints=[points[0].index, points[1].index]
	boundpoints.sort()
	to_meas_average=plot.vectors[whatset][1][boundpoints[0]:boundpoints[1]]
	avg_meas=np.mean(to_meas_average)
	std_avg_meas=np.std(to_meas_average)/np.sqrt(len(to_meas_average))
        forcebase=abs(avg_base-avg_meas)
	error=(std_avg_meas**2 + std_avg_base**2)**0.5  # pitagora error 
        print "%.1f +/- %.2f pN"%(10**12*forcebase,10**12*error)


    def plotmanip_multiplier(self, plot, current):
        '''
        Multiplies all the Y values of an SMFS curve by a value stored in the 'force_multiplier'
        configuration variable. Useful for calibrations and other stuff.
        '''
        
        #not a smfs curve...
        if current.curve.experiment != 'smfs':
            return plot
        
        #only one set is present...
        if len(self.plots[0].vectors) != 2:
            return plot
        
        #multiplier is 1...
        if (self.config['force_multiplier']==1):
            return plot

        for i in range(len(plot.vectors[0][1])):
            plot.vectors[0][1][i]=plot.vectors[0][1][i]*self.config['force_multiplier']        

        for i in range(len(plot.vectors[1][1])):
            plot.vectors[1][1][i]=plot.vectors[1][1][i]*self.config['force_multiplier']

        return plot            
   
    
    def plotmanip_flatten(self, plot, current, customvalue=False):
        '''
        Subtracts a polynomial fit to the non-contact part of the curve, as to flatten it.
        the best polynomial fit is chosen among polynomials of degree 1 to n, where n is 
        given by the configuration file or by the customvalue.
        
        customvalue= int (>0) --> starts the function even if config says no (default=False)
        '''
        
        #not a smfs curve...
        if current.curve.experiment != 'smfs':
            return plot
        
        #only one set is present...
        if len(self.plots[0].vectors) != 2:
            return plot
        
        #config is not flatten, and customvalue flag is false too
        if (not self.config['flatten']) and (not customvalue):
            return plot
        
        max_exponent=12
        delta_contact=0
        
        if customvalue:
            max_cycles=customvalue
        else:
            max_cycles=self.config['flatten'] #Using > 1 usually doesn't help and can give artefacts. However, it could be useful too.
        
        contact_index=self.find_contact_point()
        
        valn=[[] for item in range(max_exponent)]
        yrn=[0.0 for item in range(max_exponent)]
        errn=[0.0 for item in range(max_exponent)]
        
        #Check if we have a proper numerical value
        try:
            zzz=int(max_cycles)
        except:
            #Loudly and annoyingly complain if it's not a number, then fallback to zero
            print '''Warning: flatten value is not a number!
            Use "set flatten" or edit hooke.conf to set it properly
            Using zero.'''
            max_cycles=0
        
        for i in range(int(max_cycles)):
            
            x_ext=plot.vectors[0][0][contact_index+delta_contact:]
            y_ext=plot.vectors[0][1][contact_index+delta_contact:]
            x_ret=plot.vectors[1][0][contact_index+delta_contact:]
            y_ret=plot.vectors[1][1][contact_index+delta_contact:]
            for exponent in range(max_exponent):
                try:
                    valn[exponent]=sp.polyfit(x_ext,y_ext,exponent)
                    yrn[exponent]=sp.polyval(valn[exponent],x_ret)
                    errn[exponent]=sp.sqrt(sum((yrn[exponent]-y_ext)**2)/float(len(y_ext)))
                except Exception,e:
                    print 'Cannot flatten!'
                    print e
                    return plot

            best_exponent=errn.index(min(errn))
            
            #extension
            ycorr_ext=y_ext-yrn[best_exponent]+y_ext[0] #noncontact part
            yjoin_ext=np.array(plot.vectors[0][1][0:contact_index+delta_contact]) #contact part        
            #retraction
            ycorr_ret=y_ret-yrn[best_exponent]+y_ext[0] #noncontact part
            yjoin_ret=np.array(plot.vectors[1][1][0:contact_index+delta_contact]) #contact part
                
            ycorr_ext=np.concatenate((yjoin_ext, ycorr_ext))
            ycorr_ret=np.concatenate((yjoin_ret, ycorr_ret))
        
            plot.vectors[0][1]=list(ycorr_ext)
            plot.vectors[1][1]=list(ycorr_ret)
        
        return plot
            
    #---SLOPE---
    def do_slope(self,args):
        '''
        SLOPE
        (generalvclamp.py)
        Measures the slope of a delimited chunk on the return trace.
        The chunk can be delimited either by two manual clicks, or have
        a fixed width, given as an argument.
        ---------------
        Syntax: slope [width]
                The facultative [width] parameter specifies how many
                points will be considered for the fit. If [width] is
                specified, only one click will be required.
        (c) Marco Brucale, Massimo Sandal 2008
        '''

        # Reads the facultative width argument
        try:
            fitspan=int(args)
        except:
            fitspan=0

        # Decides between the two forms of user input, as per (args)
        if fitspan == 0:
            # Gets the Xs of two clicked points as indexes on the current curve vector
            print 'Click twice to delimit chunk'
            points=self._measure_N_points(N=2,whatset=1)
        else:
            print 'Click once on the leftmost point of the chunk (i.e.usually the peak)'
            points=self._measure_N_points(N=1,whatset=1)
            
        slope=self._slope(points,fitspan)

        # Outputs the relevant slope parameter
        print 'Slope:'
        print str(slope)
        to_dump='slope '+self.current.path+' '+str(slope)
        self.outlet.push(to_dump)

    def _slope(self,points,fitspan):
		# Calls the function linefit_between
        parameters=[0,0,[],[]]
        try:
            clickedpoints=[points[0].index,points[1].index]
            clickedpoints.sort()
        except:
            clickedpoints=[points[0].index-fitspan,points[0].index]        

        try:
            parameters=self.linefit_between(clickedpoints[0],clickedpoints[1])
        except:
            print 'Cannot fit. Did you click twice the same point?'
            return
             
        # Makes a vector with the fitted parameters and sends it to the GUI
        xtoplot=parameters[2]
        ytoplot=[]
        x=0
        for x in xtoplot:
            ytoplot.append((x*parameters[0])+parameters[1])
            
        clickvector_x, clickvector_y=[], []
        for item in points:
            clickvector_x.append(item.graph_coords[0])
            clickvector_y.append(item.graph_coords[1])

        lineplot=self._get_displayed_plot(0) #get topmost displayed plot
        
        lineplot.add_set(xtoplot,ytoplot)
        lineplot.add_set(clickvector_x, clickvector_y)
                
        
        if lineplot.styles==[]:
            lineplot.styles=[None,None,None,'scatter']
        else:
            lineplot.styles+=[None,'scatter']
        if lineplot.colors==[]:
            lineplot.colors=[None,None,'black',None]
        else:
            lineplot.colors+=['black',None]
        
        
        self._send_plot([lineplot])

        return parameters[0]	


    def linefit_between(self,index1,index2,whatset=1):
        '''
        Creates two vectors (xtofit,ytofit) slicing out from the
        current return trace a portion delimited by the two indexes
        given as arguments.
        Then does a least squares linear fit on that slice.
        Finally returns [0]=the slope, [1]=the intercept of the
        fitted 1st grade polynomial, and [2,3]=the actual (x,y) vectors
        used for the fit.
        (c) Marco Brucale, Massimo Sandal 2008
        '''
        # Translates the indexes into two vectors containing the x,y data to fit
        xtofit=self.plots[0].vectors[whatset][0][index1:index2]
        ytofit=self.plots[0].vectors[whatset][1][index1:index2]
        
        # Does the actual linear fitting (simple least squares with numpy.polyfit)
        linefit=[]
        linefit=np.polyfit(xtofit,ytofit,1)

        return (linefit[0],linefit[1],xtofit,ytofit)
    
    
    
