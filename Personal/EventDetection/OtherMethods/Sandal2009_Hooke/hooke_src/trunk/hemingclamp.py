#!/usr/bin/env python

'''
libhemingclamp.py

Library for interpreting Hemingway force spectroscopy files.

Copyright (C) 2008 Massimo Sandal, Marco Brucale (University of Bologna, Italy) 

This program is released under the GNU General Public License version 2.
'''
__version__='2007_02_15_devel'

__changelog__='''
2007_02_15: fixed time counter with my counter
2007_02_07: Initial implementation
'''
import string
import libhookecurve as lhc 

class DataChunk(list):
    '''Dummy class to provide ext and ret methods to the data list.
    In this case ext and self can be equal.
    '''
    
    def ext(self):
        return self
        
    def ret(self):
        return self

class hemingclampDriver(lhc.Driver):
    
    def __init__(self, filename):
        
        self.filedata = open(filename,'r')
        self.data = self.filedata.readlines()[6:]
        self.filedata.close()
        
        self.filetype = 'hemingclamp'
        self.experiment = 'clamp'
        
        self.filename=filename
       
    def __del__(self):
        self.filedata.close()   
    
    def is_me(self):
        '''
        we define our magic heuristic for HemingClamp files
        '''
        myfile=file(self.filename)
        headerlines=myfile.readlines()[0:3]
        myfile.close()
        if headerlines[0][0:10]=='#Hemingway' and headerlines[1][0:19]=='#Experiment: FClamp':
            return True
        else:
            return False
        
    def _getdata_all(self):
        time = []
        phase = []
        zpiezo = []
        defl = []
        imposed = []
        trim_indexes = []
        trim_counter = 0.0
                        
        for i in self.data:
            temp = string.split(i)
            #time.append(float(temp[0])*(1.0e-3)) # This is managed differently now, since each data point = 1ms: see below
            phase.append(float(temp[1])*(1.0e-7)) # The nonsensical (e-7) multiplier is just there to make phase data nicely plottable along other data
            zpiezo.append(float(temp[2])*(1.0e-9))
            defl.append(float(temp[3])*(1.0e-9))
            imposed.append(float(temp[4])*(1.0e-9))

        for x in range (0,len(phase)):
            if phase[x] != trim_counter:
                trim_indexes.append(x)
                trim_counter = phase[x]
       
        #we rebuild the time counter assuming 1 point = 1 millisecond
        c=0.0
        for z in zpiezo:
            time.append(c)
            c+=(1.0e-3)            
            
        return time,phase,zpiezo,defl,imposed,trim_indexes
        
    def time(self):
        return DataChunk(self._getdata_all()[0])

    def phase(self):
        return DataChunk(self._getdata_all()[1])
    
    def zpiezo(self):
        return DataChunk(self._getdata_all()[2])
     
    def deflection(self):
        return DataChunk(self._getdata_all()[3])

    def imposed(self):
        return DataChunk(self._getdata_all()[4])

    def trimindexes(self):
        return DataChunk(self._getdata_all()[5])
    
    def close_all(self):
        '''
        Explicitly closes all files
        '''
        self.filedata.close()
    
    def default_plots(self):
        main_plot=lhc.PlotObject()
        defl_plot=lhc.PlotObject()
        
        time=self.time()
        phase=self.phase()
        zpiezo=self.zpiezo()
        deflection=self.deflection()
        imposed=self.imposed()
                
        main_plot.vectors=[[time,zpiezo],[time,phase]]
        main_plot.units=['seconds','meters']
        main_plot.destination=0
        main_plot.title=self.filename
        
        defl_plot.vectors=[[time,deflection],[time,imposed]]
        defl_plot.units=['seconds','Newtons']
        defl_plot.destination=1
 
        return [main_plot, defl_plot]
    