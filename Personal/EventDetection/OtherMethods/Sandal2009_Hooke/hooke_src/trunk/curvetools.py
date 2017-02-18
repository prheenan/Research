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


class curvetoolsCommands:

      def print_prec(self, arr, prec):
         try:
           nparr=np.array(arr)
           np.set_printoptions(precision=prec)
	   #we remove the parentesis if the array is longer that 1
	   if len(nparr)>1:
  	     strvals=str(nparr)[1:-1]
           return strvals
         except:
           return "Error in the array."

      def fit_interval_nm(self,start_index,plot,nm,backwards):
	  '''
	  Calculates the number of points to fit, given a fit interval in nm
	  start_index: index of point
	  plot: plot to use
	  backwards: if true, finds a point backwards.
	  '''
	  whatset=1 #FIXME: should be decidable
	  x_vect=plot.vectors[1][0]
	  
	  c=0
	  i=start_index
	  start=x_vect[start_index]
	  maxlen=len(x_vect)
	  while abs(x_vect[i]-x_vect[start_index])*(10**9) < nm:
	      if i==0 or i==maxlen-1: #we reached boundaries of vector!
		  return c
	      
	      if backwards:
		  i-=1
	      else:
		  i+=1
	      c+=1
	  return c



      def find_current_peaks(self,noflatten, a=True, maxpeak=True, nocontact=False):
	    #Find peaks.
	    if a==True:
		  a=self.convfilt_config['mindeviation']
	    try:
		  abs_devs=float(a)
	    except:
		  print "Bad input, using default."
		  abs_devs=self.convfilt_config['mindeviation']

	    defplot=self.current.curve.default_plots()[0]
	    if not noflatten:
		flatten=self._find_plotmanip('flatten') #Extract flatten plotmanip
		defplot=flatten(defplot, self.current, customvalue=1) #Flatten curve before feeding it to has_peaks
	    pk_location,peak_size=self.has_peaks(defplot, abs_devs, maxpeak, 10, nocontact)
	    return pk_location, peak_size


      def pickup_contact_point(self,N=1,whatset=1):
	    '''macro to pick up the contact point by clicking'''
	    contact_point=self._measure_N_points(N=1, whatset=1)[0]
	    contact_point_index=contact_point.index
	    self.wlccontact_point=contact_point
	    self.wlccontact_index=contact_point.index
	    self.wlccurrent=self.current.path
	    return contact_point, contact_point_index



      def baseline_points(self,peak_location, displayed_plot):
            clicks=self.config['baseline_clicks']
            if clicks==0:
                self.basepoints=[]
                base_index_0=peak_location[-1]+self.fit_interval_nm(peak_location[-1], displayed_plot, self.config['auto_right_baseline'],False)
                self.basepoints.append(self._clickize(displayed_plot.vectors[1][0],displayed_plot.vectors[1][1],base_index_0))
                base_index_1=self.basepoints[0].index+self.fit_interval_nm(self.basepoints[0].index, displayed_plot, self.config['auto_left_baseline'],False)
                self.basepoints.append(self._clickize(displayed_plot.vectors[1][0],displayed_plot.vectors[1][1],base_index_1))
            elif clicks>0:
                print 'Select baseline'
                if clicks==1:
                    self.basepoints=self._measure_N_points(N=1, whatset=1)
                    base_index_1=self.basepoints[0].index+self.fit_interval_nm(self.basepoints[0].index, displayed_plot, self.config['auto_left_baseline'], False)
                    self.basepoints.append(self._clickize(displayed_plot.vectors[1][0],displayed_plot.vectors[1][1],base_index_1))
                else:
                    self.basepoints=self._measure_N_points(N=2, whatset=1)
            
            self.basecurrent=self.current.path
            return self.basepoints



