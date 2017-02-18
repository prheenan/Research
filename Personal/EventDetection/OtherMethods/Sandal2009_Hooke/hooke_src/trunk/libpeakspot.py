#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
a library of helping functions for spotting force spectroscopy peaks

(c)Fabrizio Benedetti and Massimo Sandal , 2007
'''
import numpy as np

def conv_dx(data,vect):
    '''
    Returns the right-centered convolution of data with vector vect
    '''
    dim=len(data)
    window=len(vect)
    temparr=[0.0]*dim
    
    end=dim-window

    for j in range(end):
     for k in range(window):
      temparr[j]+=data[j+k]*vect[k]

    return temparr  
 
def absdev(arr):
    '''
    Calculates the absolute deviation of a vector
    '''
    med=0.0
    absD=0.0
    
    med=np.mean(arr)
    for j in arr:
        absD+=abs(j-med)

    try:
      return absD/len(arr)
    except:
      return 0

 
def noise_absdev(data,positive=False,maxcut=0.2,stable=0.005):
    '''
    Returns the standard deviation of the noise.
    The logic is: we cut the most negative (or positive) data points until the absolute deviation
    becomes stable (it doesn't vary more than 0.005) or we have cut more than maxcut*len(data) points.
    Then calculate the absolute deviation.
    
    If positive=True we cut the most positive data points, False=we cut the negative ones.
    '''
    out=[item for item in data] #we copy just to be sure...
    out.sort()
    if positive:
        out.reverse()
        
    temp_absdev=absdev(out)
    cut_absdev=0
    for index in range(len(out)):
        cutindex=(index+1)*5
        cut_absdev=absdev(out[cutindex:]) #we jump five points after five points...
        if 1-(cut_absdev/temp_absdev) < stable and cutindex<(maxcut*len(out)):
            temp_absdev=cut_absdev
        else:
            break
        
    return cut_absdev
        
def abovenoise(convoluted,noise_level,cut_index=0,abs_devs=4):
    '''
    Generates a vector which is 0 where the vector is less than abs_devs*noise_level ; 1 if not (spike).
    '''
    #calculate absolute noise deviation
    #noise_level=noise_absdev(convoluted[cut_index:])
        
    above=[]
        
    for index in range(len(convoluted)):
        if index<cut_index:
            above.append(0)
        else:
            #FIXME: should calculate the *average* (here we assume that convolution mean is 0, which is FALSE!)
            if convoluted[index] < -1*noise_level*abs_devs:
                above.append(convoluted[index])
            else:
                above.append(0)
    return above
        
def find_peaks(above, seedouble=10):
    '''
    Finds individual peak location.
    abovenoise() finds all points above a given threshold in the convolution. This point is often not unique
    but there is a "cluster" of points above a threshold.
    Here we obtain only the coordinates of the largest *convolution* spike for each cluster.       
        
    above=vector obtained from abovenoise()
    seedouble=value at which we want to "delete" double peaks. That is, if two peaks have a distance
    < than $seedouble points , only the first is kept.
    '''
    nonzero=[]
    peaks_location=[]
    peaks_size=[]
        
    for index in range(len(above)):
        if above[index] != 0:
            nonzero.append(index)
        else:
            if len(nonzero) != 0:
                if len(nonzero)==1:
                    nonzero.append(nonzero[0]+1)
                peaks_size.append(min(above[nonzero[0]:nonzero[-1]]))
                peaks_location.append(above[nonzero[0]:nonzero[-1]].index(peaks_size[-1])+nonzero[0])
                nonzero=[]
            else:
                pass
                
    #recursively eliminate double peaks
    #(maybe not the smartest of cycles, but i'm asleep...)
    temp_location=None
    while temp_location != []:
        temp_location=[]
        if len(peaks_location)>1:
            for index in range(len(peaks_location)-1):
                if peaks_location[index+1]-peaks_location[index] < seedouble:
                    temp_location=peaks_location[:index]+peaks_location[index+1:]
        if temp_location != []:
            peaks_location=temp_location           
        
    return peaks_location,peaks_size