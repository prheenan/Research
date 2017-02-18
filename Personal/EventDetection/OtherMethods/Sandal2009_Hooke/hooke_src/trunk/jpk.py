#!/usr/bin/env python

import string
import libhookecurve as lhc 

class DataChunk(list):
    #Dummy class to provide ext and ret methods to the data list.
    
    def ext(self):
        halflen=(len(self)/2)
        return self[0:halflen]
        
    def ret(self):
        halflen=(len(self)/2)
        return self[halflen:]

class jpkDriver(lhc.Driver):

    def __init__(self, filename):
        self.filename=filename #self.filename can always be useful, and should be defined
        self.filedata = open(filename,'r') #We open the file
        self.filelines=self.filedata.readlines()
        self.filedata.close()        
        '''These are two strings that can be used by Hooke commands/plugins to understand what they are looking at. They have no other
        meaning. They have to be somehow defined however - commands often look for those variables.
        
        self.filetype should contain the name of the exact filetype defined by the driver (so that filetype-specific commands can know
                      if they're dealing with the correct filetype)
        self.experiment should contain instead the type of data involved (for example, various drivers can be used for force-clamp experiments,
                      but hooke commands could like to know if we're looking at force clamp data, regardless of their origin, and not other 
                      kinds of data)
        
        Of course, all other variables you like can be defined in the class.
        '''
        self.filetype = 'jpk'
        self.experiment = 'smfs'
                  
        

    def __del__(self):
        self.filedata.close()   
    
    def is_me(self):
        '''
        we define our magic heuristic for jpk files
        '''
        myfile=file(self.filename)
        headerlines=myfile.readlines()[0:3]
        myfile.close()
        if headerlines[0][0:11]=='# xPosition' and headerlines[1][0:11]=='# yPosition':
            return True
        else:
            return False
        
    def close_all(self):
        self.filedata.close()
        
    def _read_data_segment(self):
        #routine that actually reads the data
        
        height_ms=[]
        height_m=[]
        height=[]
        v_deflection=[]
        h_deflection=[]
        
        self.springconstant=0 #if we don't meet any spring constant, use deflection...
        
        for line in self.filelines:
            #we meet the segment defining the order of data columns
            
            if line[0:9]=='# columns':
                splitline=line.split()[2:]
                height_ms_index=splitline.index('smoothedStrainGaugeHeight')
                height_m_index=splitline.index('strainGaugeHeight')
                height_index=splitline.index('height')
                v_deflection_index=splitline.index('vDeflection')
                #h_deflection=splitline.index('hDeflection')
            
            if line[0:16]=='# springConstant':
                self.springconstant=float(line.split()[2])
    
            if line[0] != '#' and len(line.split())>1:
                dataline=line.split()
                height_ms.append(float(dataline[height_ms_index]))
                height_m.append(float(dataline[height_m_index]))
                height.append(float(dataline[height_index]))
                v_deflection.append(float(dataline[v_deflection_index]))
                #h_deflection.append(float(dataline[h_deflection_index]))
        
        if self.springconstant != 0:
            force=[item*self.springconstant for item in v_deflection]
        else: #we have measured no spring constant :(
            force=v_deflection
            
        height_ms=DataChunk([item*-1 for item in height_ms])
        height_m=DataChunk([item*-1 for item in height_m])
        height=DataChunk([item*-1 for item in height])
        deflection=DataChunk(v_deflection)
        force=DataChunk(force)
        
        return height_ms,height_m,height,deflection,force
        
    def deflection(self):
        height_ms,height_m,height,deflection,force=self._read_data_segment()
        deflection_ext=deflection.ext()
        deflection_ret=deflection.ret()
        deflection_ret.reverse()
        return deflection_ext,deflection_ret
        
    def default_plots(self):
        
        height_ms,height_m,height,deflection,force=self._read_data_segment()
        
        height_ms_ext=height_ms.ext()
        height_ms_ret=height_ms.ret()
        force_ext=force.ext()
        force_ret=force.ret()
        #reverse the return data, to make it coherent with hooke standard
        height_ms_ret.reverse()
        force_ret.reverse()
        
        main_plot=lhc.PlotObject()  
        main_plot.add_set(height_ms_ext,force_ext)
        main_plot.add_set(height_ms_ret,force_ret)
        
        
        
        if self.springconstant != 0:
            main_plot.units=['meters','force']
        else:
            main_plot.units=['meters','meters']
        
        main_plot.normalize_vectors()
        
        main_plot.destination=0
        main_plot.title=self.filename
        
        return [main_plot]