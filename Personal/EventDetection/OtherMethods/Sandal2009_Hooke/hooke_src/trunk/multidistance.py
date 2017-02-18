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


class multidistanceCommands:
    
    def do_multidistance(self,args):
     '''
     MULTIDISTANCE
     multidistance.py
     Based on the convolution recognition automatically give the distances between the peaks found.
     The command allow also to remove the unwanted peaks that can be due to interference.
     When you first issue the command, it will ask for the filename. If you are giving the filename
     of an existing file, autopeak will resume it and append measurements to it. If you are giving
     a new filename, it will create the file and append to it until you close Hooke.
     You can also define a minimun deviation of the peaks.
     
     Syntax:
     multidistance [deviation]
     deviation = number of times the convolution signal is above the noise absolute deviation.
     '''

      
     noflatten=False
     peaks_location, peak_size=self.find_current_peaks(noflatten)
     
     #if no peaks, we have nothing to plot. exit.
     if len(peaks_location)==0:
            return
        
     #otherwise, we plot the peak locations.
     xplotted_ret=self.plots[0].vectors[1][0]
     yplotted_ret=self.plots[0].vectors[1][1]
     xgood=[xplotted_ret[index] for index in peaks_location]
     ygood=[yplotted_ret[index] for index in peaks_location]
        
     recplot=self._get_displayed_plot()
     recplot.vectors.append([xgood,ygood])
     if recplot.styles==[]:
         recplot.styles=[None,None,'scatter']
         recplot.colors=[None,None,None]
     else:
         recplot.styles+=['scatter']
         recplot.colors+=[None]

     self._send_plot([recplot])

     print 'Peaks to ignore (0,1...n from contact point,return to take all)'
     print 'N to discard measurement'
     exclude_raw=raw_input('Input:')
     if exclude_raw=='N':
        print 'Discarded.'
        return
     
     if not exclude_raw=='':
        exclude=exclude_raw.split(',')
	#we convert in numbers the input
        try:
            exclude=[int(item) for item in exclude]
        except:
            print 'Bad input, taking nothing.'
	    return

#    we remove the peaks that we don't want from the list, we need a counter beacause if we remove
#    a peaks the other peaks in the list are shifted by one at each step
        count=0
        for a in exclude:
          if (a==0):
             peaks_location=peaks_location[1:]
          else:
             new_a=a-count
             peaks_location=  peaks_location[0:new_a]+peaks_location[new_a+1:]
             peak_size=            peak_size[0:new_a]+peak_size[new_a+1:]
          count+=1
     
     #we calculate the distance vector
     dist=[]
     for i in range(len(peaks_location)-1):
         dist.append(xplotted_ret[peaks_location[i]]-xplotted_ret[peaks_location[i+1]])
     
     



        #Save file info
     if self.autofile=='':
            self.autofile=raw_input('Multidistance filename? (return to ignore) ')
            if self.autofile=='':
                print 'Not saved.'
                return
        
     if not os.path.exists(self.autofile):
            f=open(self.autofile,'w+')
            f.write('Analysis started '+time.asctime()+'\n')
            f.write('----------------------------------------\n')
            f.write('; Delta Distance length (m)\n')
	    f.write(self.current.path+'\n')
	    for o in dist:
	       f.write(";")
	       f.write(str(o))
            f.write("\n")
            f.close()
            
     print 'Saving...'
     f=open(self.autofile,'a+')
        
     f.write(self.current.path+'\n')
     for i in dist:
          f.write(";")
          f.write(str(i))

     f.write("\n")            
     f.close()
