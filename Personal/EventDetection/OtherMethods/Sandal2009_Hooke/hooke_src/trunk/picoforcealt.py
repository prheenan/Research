#!/usr/bin/env python

'''
libpicoforce.py

Library for interpreting Picoforce force spectroscopy files. Alternate version

Copyright (C) 2006 Massimo Sandal (University of Bologna, Italy).
Copyright (C) 2008 Alberto Gomez-Casado (University of Twente, Netherlands).

This program is released under the GNU General Public License version 2.
'''

import re, struct
from scipy import arange

import libhookecurve as lhc

__version__='0.0.0.20081706'



class DataChunk(list):
    #Dummy class to provide ext and ret methods to the data list.
    
    def ext(self):
        halflen=(len(self)/2)
        return self[0:halflen]
        
    def ret(self):
        halflen=(len(self)/2)
        return self[halflen:]

class picoforcealtDriver(lhc.Driver):

    #Construction and other special methods
    
    def __init__(self,filename):
        '''
        constructor method
        '''
        
        self.textfile=file(filename)
        self.binfile=file(filename,'rb')
        
        #The 0,1,2 data chunks are:
        #0: D (vs T)
        #1: Z (vs T)
        #2: D (vs Z)
        
        self.forcechunk=0
        self.distancechunk=1
	#TODO eliminate the need to set chunk numbers
        
        self.filepath=filename
        self.debug=True
        
        self.filetype='picoforce'
        self.experiment='smfs'
    
    
            
    def _get_samples_line(self):
        '''
        Gets the samples per line parameters in the file, to understand trigger behaviour.
        '''
        self.textfile.seek(0)
        
        samps_expr=re.compile(".*Samps")
        
        samps_values=[]
        for line in self.textfile.readlines():
            if samps_expr.match(line):
                try:
                    samps=int(line.split()[2]) #the third word splitted is the offset (in bytes)
                    samps_values.append(samps)
                except:
                    pass
                
                #We raise a flag for the fact we meet an offset, otherwise we would take spurious data length arguments.
                        
        return int(samps_values[0])
                    
    def _get_chunk_coordinates(self):
        '''
        This method gets the coordinates (offset and length) of a data chunk in our
        Picoforce file.
        
        It returns a list containing two tuples: 
        the first element of each tuple is the data_offset, the second is the corresponding 
        data size.
        
        In near future probably each chunk will get its own data structure, with 
        offset, size, type, etc.
        '''
        self.textfile.seek(0)
        
        offset_expr=re.compile(".*Data offset")
        length_expr=re.compile(".*Data length")

        data_offsets=[]
        data_sizes=[]
        flag_offset=0

        for line in self.textfile.readlines():

            if offset_expr.match(line):
                offset=int(line.split()[2]) #the third word splitted is the offset (in bytes)
                data_offsets.append(offset)
                #We raise a flag for the fact we meet an offset, otherwise we would take spurious data length arguments.
                flag_offset=1 
    
            #same for the data length
            if length_expr.match(line) and flag_offset: 
                size=int(line.split()[2])
                data_sizes.append(size)
                #Put down the offset flag until the next offset is met.
                flag_offset=0

        return zip(data_offsets,data_sizes)
        
    def _get_data_chunk(self,whichchunk):
        '''
        reads a data chunk and converts it in 16bit signed int.
        '''
        offset,size=self._get_chunk_coordinates()[whichchunk]
        
        
        self.binfile.seek(offset)
        raw_chunk=self.binfile.read(size)
        
        my_chunk=[]
        for data_position in range(0,len(raw_chunk),2):
            data_unit_bytes=raw_chunk[data_position:data_position+2]
            #The unpack function converts 2-bytes in a signed int ('h').
            #we use output[0] because unpack returns a 1-value tuple, and we want the number only
            data_unit=struct.unpack('h',data_unit_bytes)[0]
            my_chunk.append(data_unit)                             
        
        return DataChunk(my_chunk)

    def _force(self):
	#returns force vector
        Kspring=self.get_spring_constant()
        return DataChunk([(meter*Kspring) for meter in self._deflection()])

    def _deflection(self):
        #for internal use (feeds _force)
        voltrange=1
        z_scale=self._get_Z_scale()
        deflsensitivity=self.get_deflection_sensitivity()
        volts=[((float(lsb))*voltrange*z_scale) for lsb in self.data_chunks[self.forcechunk]]    
        deflect=[volt*deflsensitivity for volt in volts]
   	     
        return deflect
             
    
    def _Z(self):
        #returns distance vector (calculated instead than from data chunk)        
        rampsize=self._get_rampsize()
        sampsline=self._get_samples_line()
        senszscan=self._get_Z_scan_sens()
             
        xstep=senszscan*rampsize/sampsline*10**(-9)

        xext=arange(sampsline*xstep,0,-xstep)
        xret=arange(sampsline*xstep,0,-xstep)
         
        return DataChunk(xext.tolist()+xret.tolist())
    
    def _get_Z_scale(self):
        self.textfile.seek(0)
        expr=re.compile(".*@4:Z scale")
        
        for line in self.textfile.readlines():
            if expr.match(line):
                zscale=float((line.split()[5]).strip("() []"))
                break
        return zscale
    
    def _get_rampsize(self):
        self.textfile.seek(0)
        expr=re.compile(".*@4:Ramp size:")
        
        for line in self.textfile.readlines():
            if expr.match(line):
                zsens=float((line.split()[7]).strip("() []"))
                break
        return zsens
        
    def _get_Z_scan_sens(self):
        self.textfile.seek(0)
        expr=re.compile(".*@Sens. Zsens")
        
        for line in self.textfile.readlines():
            if expr.match(line):
                zsens=float((line.split()[3]).strip("() []"))
                break
        return zsens
    
    
                
    def get_deflection_sensitivity(self):
        '''
        gets deflection sensitivity
        '''    
        self.textfile.seek(0)
        
        def_sensitivity_expr=re.compile(".*@Sens. DeflSens")
        
        for line in self.textfile.readlines():
            if def_sensitivity_expr.match(line):
                def_sensitivity=float(line.split()[3])
                break
        #return it in SI units (that is: m/V, not nm/V)
        return def_sensitivity*(10**(-9))
        
    def get_spring_constant(self):
        '''
        gets spring constant.
        We actually find *three* spring constant values, one for each data chunk (F/t, Z/t, F/z).
        They are normally all equal, but we retain all three for future...
        '''
        self.textfile.seek(0)
        
        springconstant_expr=re.compile(".*Spring Constant")
        
        constants=[]
        
        for line in self.textfile.readlines():
            if springconstant_expr.match(line):
                constants.append(float(line.split()[2]))
        
        return constants[0]
    
    def is_me(self):
        '''
        self-identification of file type magic
        '''
        curve_file=file(self.filepath)
        header=curve_file.read(30)
        curve_file.close()
        
        if header[2:17] == 'Force file list': #header of a picoforce file
            #here DONT translate chunk
            self.data_chunks=[self._get_data_chunk(num) for num in [0,1,2]]
            return True
        else:
            return False
    
    def close_all(self):
        '''
        Explicitly closes all files
        '''
        self.textfile.close()
        self.binfile.close()
    
    def default_plots(self):
        '''
        creates the default PlotObject
        '''
        force=self._force()
        zdomain=self._Z()
        samples=self._get_samples_line()
        main_plot=lhc.PlotObject()
        main_plot.vectors=[[zdomain.ext()[0:samples], force.ext()[0:samples]],[zdomain.ret()[0:samples], force.ret()[0:samples]]]
        main_plot.normalize_vectors()
        main_plot.units=['meters','newton']
        main_plot.destination=0
        main_plot.title=self.filepath
        
        
        return [main_plot]

    def deflection(self):
        #interface for correct plotmanip and others
        deflectionchunk=DataChunk(self._deflection())
        return deflectionchunk.ext(),deflectionchunk.ret()
