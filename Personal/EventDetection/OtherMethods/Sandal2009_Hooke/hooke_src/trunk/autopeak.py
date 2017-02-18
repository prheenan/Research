#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


class autopeakCommands:
    
    def do_autopeak(self,args):
        '''
        AUTOPEAK
        (autopeak.py)
        Automatically performs a number of analyses on the peaks of the given curve.
        Currently it automatically:
        - fits peaks with WLC or FJC function (depending on how the fit_function variable is set)
        - measures peak maximum forces with a baseline
        - measures slope in proximity of peak maximum
        Requires flatten plotmanipulator , fit.py plugin , flatfilts.py plugin with convfilt
        
        Syntax:
        autopeak [rebase] [pl=value] [manual=value] [t=value] [noauto] [reclick]
        
        rebase : Re-asks baseline interval
        
        pl=[value] : Use a fixed persistent length for the fit. If pl is not given, 
                     the fit will be a 2-variable  
                     fit. DO NOT put spaces between 'pl', '=' and the value.
                     The value must be in meters. 
                     Scientific notation like 0.35e-9 is fine.

        manual=[value]:  Allow to choose the peaks to analyze. It need, as a value, the number
                         of the peaks to analyze. NOTE: It have to be used with the manual selection of the baseline.

        t=[value] : Use a user-defined temperature. The value must be in
                    kelvins; by default it is 293 K.
                    DO NOT put spaces between 't', '=' and the value.
        
        noauto : allows for clicking the contact point by 
                 hand (otherwise it is automatically estimated) the first time.
                 If subsequent measurements are made, the same contact point
                 clicked the first time is used
        
        reclick : redefines by hand the contact point, if noauto has been used before
                  but the user is unsatisfied of the previously choosen contact point.
        
        usepoints : fit interval by number of points instead than by nanometers
        
        noflatten : does not use the "flatten" plot manipulator

	nocontact : when there is a curve without the contact region we have to put this flag

        When you first issue the command, it will ask for the filename. If you are giving the filename
        of an existing file, autopeak will resume it and append measurements to it. If you are giving
        a new filename, it will create the file and append to it until you close Hooke.
        
        
        Useful variables (to set with SET command):
        ---
        fit_function = type of function to use for elasticity. If "wlc" worm-like chain is used, if "fjc" freely jointed
                       chain is used
        
        temperature= temperature of the system for wlc/fjc fit (in K)
        
        auto_slope_span = number of points on which measure the slope, for slope
        
        auto_fit_nm = number of nm to fit before the peak maximum, for WLC/FJC (if usepoints false)
        auto_fit_points = number of points to fit before the peak maximum, for WLC/FJC (if usepoints true)
        
        baseline_clicks = -1: no baseline, f=0 at the contact point (whether hand-picked or automatically found)
                           0: automatic baseline
                           1: decide baseline with a single click and length defined in auto_left_baseline
                           2: let user click points of baseline
        auto_left_baseline = length in nm to use as baseline from the right point (if baseline_clicks=0 , 1)
        auto_right_baseline = distance in nm of peak-most baseline point from last peak (if baseline_clicks = 0)
        
        auto_min_p ; auto_max_p = Minimum and maximum persistence length (if using WLC) or Kuhn length (if using FJC)
                                outside of which the peak is automatically discarded (in nm)
        '''
            
        #default fit etc. variables
        pl_value=None
        T=self.config['temperature']
        
        slope_span=int(self.config['auto_slope_span'])
        delta_force=10
        rebase=False #if true=we select rebase
        noflatten=False #if true=we avoid flattening
	nocontact=False

        #initialize output data vectors
        c_lengths=[]
        p_lengths=[]
        sigma_c_lengths=[]
        sigma_p_lengths=[]
        forces=[]
        slopes=[]
        qstds=[]
        manualpoints=0

        #pick up plot
        displayed_plot=self._get_displayed_plot(0)
        
        #COMMAND LINE PARSING
        #--Using points instead of nm interval
        if 'usepoints' in args.split():
            fit_points=int(self.config['auto_fit_points'])
            usepoints=True
        else:
            fit_points=None
            usepoints=False
        #--Recalculate baseline
        if 'rebase' in args or (self.basecurrent != self.current.path):
            rebase=True 
        
        if 'noflatten' in args:
            noflatten=True

        if 'nocontact' in args:
            nocontact=True 
	    #print "nocontact TRUE"

        #--Custom persistent length / custom temperature
        for arg in args.split():
            #look for a persistent length argument.
            if 'pl=' in arg:
                pl_expression=arg.split('=')
                pl_value=float(pl_expression[1]) #actual value
            #look for a T argument. FIXME: spaces are not allowed between 'pl' and value
            if ('t=' in arg[0:3]) or ('T=' in arg[0:3]):
                t_expression=arg.split('=')
                T=float(t_expression[1])
            if ('manual=' in arg):
		print "Manual selection of the interesting peaks."
		manualpoints_expression=arg.split('=')
		manualpoints=int(manualpoints_expression[1])
        #--Contact point arguments
        if 'reclick' in args.split():
            print 'Click contact point'
            contact_point, contact_point_index = self.pickup_contact_point()
        elif 'noauto' in args.split():
            if self.wlccontact_index==None or self.wlccurrent != self.current.path:
                print 'Click contact point'
                contact_point , contact_point_index = self.pickup_contact_point()
            else:
                contact_point=self.wlccontact_point
                contact_point_index=self.wlccontact_index
        else:
            #Automatically find contact point
            cindex=self.find_contact_point()
            contact_point=self._clickize(displayed_plot.vectors[1][0], displayed_plot.vectors[1][1], cindex)
        #--END COMMAND LINE PARSING--
        
        

        if not manualpoints:
           peak_location, peak_size = self.find_current_peaks(noflatten, True, True, nocontact)
	else:
	   peak_location=[]
           for i in range(manualpoints):
	      print "Select manually the peaks:"
	      peak_location.append(  (self._measure_N_points(N=1, whatset=1)[0]).index  )
        
        if len(peak_location) == 0:
            print 'No peaks to fit.'
            return
        
        fitplot=copy.deepcopy(displayed_plot)
        
        #Pick up force baseline
        if rebase:
            self.basepoints=self.baseline_points(peak_location, displayed_plot)
        
        boundaries=[self.basepoints[0].index, self.basepoints[1].index]
        boundaries.sort()
        to_average=displayed_plot.vectors[1][1][boundaries[0]:boundaries[1]] #y points to average
        avg=np.mean(to_average)
        
        clicks=self.config['baseline_clicks']
        if clicks==-1:
            try:
                avg=displayed_plot.vectors[1][1][contact_point_index]
            except:
                avg=displayed_plot.vectors[1][1][cindex]
        
        for peak in peak_location:
            #WLC FITTING
            #define fit interval
            if not usepoints:
                fit_points=self.fit_interval_nm(peak, displayed_plot, self.config['auto_fit_nm'], True)
            peak_point=self._clickize(displayed_plot.vectors[1][0],displayed_plot.vectors[1][1],peak)
            other_fit_point=self._clickize(displayed_plot.vectors[1][0],displayed_plot.vectors[1][1],peak-fit_points)
            
            #points for the fit
            points=[contact_point, peak_point, other_fit_point]
            
            if abs(peak_point.index-other_fit_point.index) < 2:
                continue
            
            if self.config['fit_function']=='wlc':
                
                params, yfit, xfit, fit_errors, qstd = self.wlc_fit(points, displayed_plot.vectors[1][0], displayed_plot.vectors[1][1], pl_value, T, return_errors=True)
            elif self.config['fit_function']=='fjc':
                params, yfit, xfit, fit_errors, qstd = self.fjc_fit(points, displayed_plot.vectors[1][0], displayed_plot.vectors[1][1], pl_value, T, return_errors=True)
            else:
                print 'Unknown fit function'
                print 'Please set fit_function as wlc or fjc'
                return
            
            
            #Measure forces
            delta_to_measure=displayed_plot.vectors[1][1][peak-delta_force:peak+delta_force]
            y=min(delta_to_measure)
            #save force values (pN)   
            #Measure slopes
            slope=self.linefit_between(peak-slope_span,peak)[0]
            
            
            #check fitted data and, if right, add peak to the measurement
            #FIXME: code duplication
            if len(params)==1: #if we did choose 1-value fit
                p_lengths.append(pl_value)
                c_lengths.append(params[0]*(1.0e+9))
                sigma_p_lengths.append(0)
                sigma_c_lengths.append(fit_errors[0]*(1.0e+9))
                forces.append(abs(y-avg)*(1.0e+12))
                slopes.append(slope)
                qstds.append(qstd)     
                #Add WLC fit lines to plot
                fitplot.add_set(xfit,yfit)
                if len(fitplot.styles)==0:
                    fitplot.styles=[]
                    fitplot.colors=[]
                else:
                    fitplot.styles.append(None)
                    fitplot.colors.append(None)
            else: #2-value fit
                p_leng=params[1]*(1.0e+9)
                #check if persistent length makes sense. otherwise, discard peak.
                if p_leng>self.config['auto_min_p'] and p_leng<self.config['auto_max_p']:
                    p_lengths.append(p_leng)       
                    c_lengths.append(params[0]*(1.0e+9))
                    sigma_c_lengths.append(fit_errors[0]*(1.0e+9))
                    sigma_p_lengths.append(fit_errors[1]*(1.0e+9))
                    forces.append(abs(y-avg)*(1.0e+12))
                    slopes.append(slope)
                    qstds.append(qstd)     
                    
                    #Add WLC fit lines to plot
                    fitplot.add_set(xfit,yfit)
                    if len(fitplot.styles)==0:
                        fitplot.styles=[]
                        fitplot.colors=[]
                    else:
                        fitplot.styles.append(None)
                        fitplot.colors.append(None)
                else:
                    pass
 
            
        #add basepoints to fitplot
        fitplot.add_set([self.basepoints[0].graph_coords[0],self.basepoints[1].graph_coords[0]],[self.basepoints[0].graph_coords[1],self.basepoints[1].graph_coords[1]]) 
        fitplot.styles.append('scatter')
        fitplot.colors.append(None)
        
        #Show wlc fits and peak locations
        self._send_plot([fitplot])
        
        print 'Using fit function: ',self.config['fit_function']
        print 'Measurements for all peaks detected:'
        print 'contour (nm)', self.print_prec(c_lengths,1)
        print 'sigma contour (nm)',self.print_prec(sigma_c_lengths,2)
        print 'p (nm)',self.print_prec(p_lengths,3)
        print 'sigma p (nm)',self.print_prec(sigma_p_lengths,3)
        print 'forces (pN)',self.print_prec(forces,1)
        print 'slopes (N/m)',self.print_prec(slopes,3)
        
        controller=False
        while controller==False:
	  #Ask the user what peaks to ignore from analysis.
	  print 'Peaks to ignore (0,1...n from contact point,return to take all)'
	  print 'N to discard measurement'
	  exclude_raw=raw_input('Input:')
	  if exclude_raw=='N':
	      print 'Discarded.'
	      return
	  if exclude_raw=='':
	      controller=True
	  if not exclude_raw=='':
	      exclude=exclude_raw.split(',')
	      try:
		  exclude=[int(item) for item in exclude]
		  for i in exclude:
		      c_lengths[i]=None
		      p_lengths[i]=None
		      forces[i]=None
		      slopes[i]=None
		      sigma_c_lengths[i]=None
		      sigma_p_lengths[i]=None
                      controller=True
	      except:
		  print 'Bad input.'
                  controller=False

        #Clean data vectors from ignored peaks        
        #FIXME:code duplication
        c_lengths=[item for item in c_lengths if item != None]
        p_lengths=[item for item in p_lengths if item != None]
        forces=[item for item in forces if item != None]
        slopes=[item for item in slopes if item != None]    
        sigma_c_lengths=[item for item in sigma_c_lengths if item != None]    
        sigma_p_lengths=[item for item in sigma_p_lengths if item != None]    
        
        print 'Measurements for chosen peaks:'
        print 'contour (nm)', self.print_prec(c_lengths,1)
        print 'sigma contour (nm)',self.print_prec(sigma_c_lengths,2)
        print 'p (nm)',self.print_prec(p_lengths,3)
        print 'sigma p (nm)',self.print_prec(sigma_p_lengths,3)
        print 'forces (pN)',self.print_prec(forces,1)
        print 'slopes (N/m)',self.print_prec(slopes,3)
        
        #Save file info
        if self.autofile=='':
            self.autofile=raw_input('Autopeak filename? (return to ignore) ')
            if self.autofile=='':
                print 'Not saved.'
                return
        
        if not os.path.exists(self.autofile):
            f=open(self.autofile,'w+')
            f.write('Analysis started '+time.asctime()+'\n')
            f.write('----------------------------------------\n')
            f.write('; Contour length (nm)  ;  Persistence length (nm) ;  Max.Force (pN)  ;  Slope (N/m) ;  Sigma contour (nm) ; Sigma persistence (nm)\n')
            f.close()
            
        print 'Saving...'
        f=open(self.autofile,'a+')
        
        f.write(self.current.path+'\n')
        for i in range(len(c_lengths)):
            f.write(' ; '+str(c_lengths[i])+' ; '+str(p_lengths[i])+' ; '+str(forces[i])+' ; '+str(slopes[i])+' ; '+str(sigma_c_lengths[i])+' ; '+str(sigma_p_lengths[i])+'\n')
            
        f.close()
        self.do_note('autopeak')
        
