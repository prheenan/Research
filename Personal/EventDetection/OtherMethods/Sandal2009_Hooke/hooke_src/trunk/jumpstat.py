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
import sys
import warnings
warnings.simplefilter('ignore',np.RankWarning)


class jumpstatCommands():
    
    def do_jumpstat(self,args):
     '''
     JUMPSTAT
     jumpstat.py
     Based on the convolution recognition automatically give:
     - the delta distance between the peaks,
     - the delta-force from the top of the peaks and subsequent relaxation,
     - the delta-force from the top of the peaks and the baseline
     The command allow also to remove the unwanted peaks that can be due to interference.
     When you first issue the command, it will ask for the filename. If you are giving the filename
     of an existing file, autopeak will resume it and append measurements to it. If you are giving
     a new filename, it will create the file and append to it until you close Hooke.
     You can also define a minimun deviation of the peaks.
     
     Syntax:
     jumpstat [deviation]
     deviation = number of times the convolution signal is above the noise absolute deviation.
     '''


     #finding the max and the minimum positions for all the peaks
     noflatten=False
     #we use if else only to avoid a "bad input" message from find_current_peaks
     if (len(args)==0):
             max_peaks_location, peak_size=self.find_current_peaks(noflatten)
	     min_peaks_location, pks2=self.find_current_peaks(noflatten, True, False)
     else:
	     max_peaks_location, peak_size=self.find_current_peaks(noflatten, args)
	     min_peaks_location, pks2=self.find_current_peaks(noflatten, args, False)


     #print "max_peaks_location:  "+str(len(max_peaks_location))
     #print "min_peaks_location:  "+str(len(min_peaks_location))

     #if no peaks, we have nothing to plot. exit.
     if len(max_peaks_location)==0:
            print "No peaks on this curve."
            return

     if len(max_peaks_location)!=len(min_peaks_location):
            print "Something went wrong in peaks recognition, number of minima is different from number of maxima. Exiting."
            return

     #otherwise, we plot the peak locations.
     xplotted_ret=self.plots[0].vectors[1][0]
     yplotted_ret=self.plots[0].vectors[1][1]
     xgood=[xplotted_ret[index] for index in max_peaks_location]
     ygood=[yplotted_ret[index] for index in max_peaks_location]

     xafter=[xplotted_ret[index] for index in min_peaks_location]
     yafter=[yplotted_ret[index] for index in min_peaks_location]

     recplot=self._get_displayed_plot()
     recplot2=self._get_displayed_plot()
     recplot.vectors.append([xgood,ygood])
     recplot2.vectors.append([xafter,yafter])

     if recplot.styles==[]:
         recplot.styles=[None,None,'scatter']
         recplot.colors=[None,None,None]
     else:
         recplot.styles+=['scatter']
         recplot.colors+=[None]

     if recplot2.styles==[]:
         recplot2.styles=[None,None,None]
         recplot2.colors=[None,'1.0',None]
     else:
         recplot2.styles+=['scatter']
         recplot2.colors+=['0.5']

     self._send_plot([recplot])
     self._send_plot([recplot2])


     #finding the baseline
     self.basepoints=self.baseline_points(max_peaks_location, recplot)    
     boundaries=[self.basepoints[0].index, self.basepoints[1].index]
     boundaries.sort()
     to_average=recplot.vectors[1][1][boundaries[0]:boundaries[1]] #y points to average
     avg=np.mean(to_average)


     dist=[]
     jumpforce=[]
     force=[]

     #we calculate the distance vector
     for g in range(len(max_peaks_location)-1):
       dist.append((10**9)*(xplotted_ret[max_peaks_location[g]]-xplotted_ret[max_peaks_location[g+1]]))
     print "Distance values for the peaks in nm:"
     print dist

     #the jump-force vector
     for g in range(len(max_peaks_location)):
      jumpforce.append((10**12)     *(yplotted_ret[min_peaks_location[g]] -yplotted_ret[max_peaks_location[g]])   )
     print "Force values for the jumps of the peaks in pN:"
     print jumpforce

     #the force from baseline vector
     for g in range(len(max_peaks_location)):
      force.append((10**12)*(avg-yplotted_ret[max_peaks_location[g]]))
     print "Force values for the peaks in pN:"
     print force
     


     #Now ask for the peaks that we don't want
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
             max_peaks_location=max_peaks_location[1:]
	     min_peaks_location=min_peaks_location[1:]
          else:
             new_a=a-count
             max_peaks_location=  max_peaks_location[0:new_a]+max_peaks_location[new_a+1:]
             min_peaks_location=  min_peaks_location[0:new_a]+min_peaks_location[new_a+1:]
             peak_size=            peak_size[0:new_a]+peak_size[new_a+1:]
          count+=1


     #print "max_peaks_location:  "+str(len(max_peaks_location))
     #print "min_peaks_location:  "+str(len(min_peaks_location))


     dist=[]
     jumpforce=[]
     force=[]
     #we recalculate the distances and the forces after the removing of the unwanted peaks
     for g in range(len(max_peaks_location)-1):
         dist.append(xplotted_ret[max_peaks_location[g]]-xplotted_ret[max_peaks_location[g+1]])
     for g in range(len(max_peaks_location)):
         jumpforce.append( yplotted_ret[min_peaks_location[g]] - yplotted_ret[max_peaks_location[g]]   )
     for g in range(len(max_peaks_location)):
         force.append(avg  -  yplotted_ret[max_peaks_location[g]])





        #Save file info
     if self.autofile=='':
            self.autofile=raw_input('Jumpstat filename? (return to ignore) ')
            if self.autofile=='':
                print 'Not saved.'
                return

     if not os.path.exists(self.autofile):
            f=open(self.autofile,'w+')
            f.write('Analysis started '+time.asctime()+'\n')
            f.write('----------------------------------------\n')
            f.write('; Delta Distance length (m);  Jump Force pN;  Standard Force pN\n')
	    f.write(self.current.path+'\n')
	    for k in range(len(dist)):
	       f.write(";")
	       f.write(str(dist[k])+";"+str(jumpforce[k])+";"+str(force[k])+"\n"   )
            f.write("\n")
            f.close()
            
     else:
	    f=open(self.autofile,'a+')
	    f.write(self.current.path+'\n')
	    for k in range(len(dist)):
	      f.write(";")
	      f.write(str(dist[k])+";"+str(jumpforce[k])+";"+str(force[k])+"\n"   )
	    f.write("\n")
	    f.close()
 
     print 'Saving...'