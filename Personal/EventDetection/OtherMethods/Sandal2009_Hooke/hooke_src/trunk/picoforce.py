#!/usr/bin/env python

'''
libpicoforce.py

Library for interpreting Picoforce force spectroscopy files.

Copyright (C) 2006 Massimo Sandal (University of Bologna, Italy).

This program is released under the GNU General Public License version 2.
'''

import re, struct
from scipy import arange

import libhookecurve as lhc

__version__='0.0.0.20080404'


class DataChunk(list):
    #Dummy class to provide ext and ret methods to the data list.
    
    def ext(self):
        halflen=(len(self)/2)
        return self[0:halflen]
        
    def ret(self):
        halflen=(len(self)/2)
        return self[halflen:]

class picoforceDriver(lhc.Driver):

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
        
        
        self.filepath=filename
        self.debug=False
        
        self.filetype='picoforce'
        self.experiment='smfs'
    
    
    #Hidden methods. These are meant to be used only by API functions. If needed, however,
    #they can be called just like API methods.
                    
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
             
    def _get_Zscan_info(self,index):
        '''
        gets the Z scan informations needed to interpret the data chunk.
        These info come from the general section, BEFORE individual chunk headers.
        
        By itself, the function will parse for three parameters.
        (index) that tells the function what to return when called by
        exposed API methods.
        index=0 : returns Zscan_V_LSB
        index=1 : returns Zscan_V_start
        index=2 : returns Zscan_V_size
        '''
        self.textfile.seek(0)
        
        ciaoforcelist_expr=re.compile(".*Ciao force")
        zscanstart_expr=re.compile(".*@Z scan start")
        zscansize_expr=re.compile(".*@Z scan size")
                
        ciaoforce_flag=0
        theline=0
        for line in self.textfile.readlines():
            if ciaoforcelist_expr.match(line):
                ciaoforce_flag=1 #raise a flag: zscanstart and zscansize params to read are later
 
            if ciaoforce_flag and zscanstart_expr.match(line):
                raw_Zscanstart_line=line.split()
            
            if ciaoforce_flag and zscansize_expr.match(line):
                raw_Zscansize_line=line.split()

        Zscanstart_line=[]
        Zscansize_line=[]
        for itemscanstart,itemscansize in zip(raw_Zscanstart_line,raw_Zscansize_line):
            Zscanstart_line.append(itemscanstart.strip('[]()'))
            Zscansize_line.append(itemscansize.strip('[]()'))
        
        Zscan_V_LSB=float(Zscanstart_line[6])
        Zscan_V_start=float(Zscanstart_line[8])
        Zscan_V_size=float(Zscansize_line[8])
        
        return (Zscan_V_LSB,Zscan_V_start,Zscan_V_size)[index]
    
    def _get_Z_magnify_scale(self,whichchunk):
        '''
        gets Z scale and Z magnify
        Here we get Z scale/magnify from the 'whichchunk' only.
        whichchunk=1,2,3
        TODO: make it coherent with data_chunks syntaxis (0,1,2)
        
        In future, should we divide the *file* itself into chunk descriptions and gain
        true chunk data structures?
        '''
        self.textfile.seek(0)
        
        z_scale_expr=re.compile(".*@4:Z scale")
        z_magnify_expr=re.compile(".*@Z magnify")
        
        ramp_size_expr=re.compile(".*@4:Ramp size")
        ramp_offset_expr=re.compile(".*@4:Ramp offset")
                
        occurrences=0
        found_right=0
        
        
        for line in self.textfile.readlines():
            if z_magnify_expr.match(line):
                occurrences+=1
                if occurrences==whichchunk:
                    found_right=1
                    raw_z_magnify_expression=line.split()
                else:
                    found_right=0
                    
            if found_right and z_scale_expr.match(line):
                raw_z_scale_expression=line.split()
            if found_right and ramp_size_expr.match(line):
                raw_ramp_size_expression=line.split()
            if found_right and ramp_offset_expr.match(line):
                raw_ramp_offset_expression=line.split()
                
        return float(raw_z_magnify_expression[5]),float(raw_z_scale_expression[7]), float(raw_ramp_size_expression[7]), float(raw_ramp_offset_expression[7]), float(raw_z_scale_expression[5][1:])
       
    
    #Exposed APIs.
    #These are the methods that are meant to be called from external apps.          
    
    def LSB_to_volt(self,chunknum,voltrange=20):
        '''
        Converts the LSB data of a given chunk (chunknum=0,1,2) in volts.
        First step to get the deflection and the force.
        
        SYNTAXIS:
        item.LSB_to_volt(chunknum, [voltrange])
                        
        The voltrange is by default set to 20 V.
        '''
        return DataChunk([((float(lsb)/65535)*voltrange) for lsb in self.data_chunks[chunknum]])
        
    def LSB_to_deflection(self,chunknum,deflsensitivity=None,voltrange=20):
        '''
        Converts the LSB data in deflection (meters).
        
        SYNTAXIS:
        item.LSB_to_deflection(chunknum, [deflection sensitivity], [voltrange])
        
        chunknum is the chunk you want to parse (0,1,2)
        
        The deflection sensitivity by default is the one parsed from the file. 
        The voltrange is by default set to 20 V.
        '''
        if deflsensitivity is None:
            deflsensitivity=self.get_deflection_sensitivity()
            
        lsbvolt=self.LSB_to_volt(chunknum)     
        return DataChunk([volt*deflsensitivity for volt in lsbvolt])
        
    def deflection(self):
        '''
        Get the actual force curve deflection.
        '''
        deflchunk= self.LSB_to_deflection(2)
        return deflchunk.ext(),deflchunk.ret()
        
    def LSB_to_force(self,chunknum=2,Kspring=None,voltrange=20):
        '''
        Converts the LSB data (of deflection) in force (newtons).
        
        SYNTAXIS:
        item.LSB_to_force([chunknum], [spring constant], [voltrange])
        
        chunknum is the chunk you want to parse (0,1,2). The chunk used is by default 2.
        The spring constant by default is the one parsed from the file.
        The voltrange is by default set to 20 V.
        '''
        if Kspring is None:
            Kspring=self.get_spring_constant()
            
        lsbdefl=self.LSB_to_deflection(chunknum)        
        return DataChunk([(meter*Kspring) for meter in lsbdefl])
        
    def get_Zscan_V_start(self):
        return self._get_Zscan_info(1)
        
    def get_Zscan_V_size(self):
        return self._get_Zscan_info(2)
        
    def get_Z_scan_sensitivity(self):
        '''
        gets Z sensitivity
        '''
        self.textfile.seek(0)
        
        z_sensitivity_expr=re.compile(".*@Sens. Zsens")
        
        for line in self.textfile.readlines():
            if z_sensitivity_expr.match(line):
                z_sensitivity=float(line.split()[3])
        #return it in SI units (that is: m/V, not nm/V)
        return z_sensitivity*(10**(-9))
          
    def get_Z_magnify(self,whichchunk):
        '''
        Gets the Z magnify factor. Normally it is 1, unknown exact use as of 2006-01-13
        '''
        return self._get_Z_magnify_scale(whichchunk)[0]
    
    def get_Z_scale(self,whichchunk):
        '''
        Gets the Z scale.
        '''
        return self._get_Z_magnify_scale(whichchunk)[1]
    
    def get_ramp_size(self,whichchunk):
        '''
        Gets the -user defined- ramp size
        '''
        return self._get_Z_magnify_scale(whichchunk)[2]
        
    def get_ramp_offset(self,whichchunk):
        '''
        Gets the ramp offset
        '''
        return self._get_Z_magnify_scale(whichchunk)[3]
        
    def get_Z_scale_LSB(self,whichchunk):
        '''
        Gets the LSB-to-volt conversion factor of the Z data.
        (so called hard-scale in the Nanoscope documentation)
        
        '''
        return self._get_Z_magnify_scale(whichchunk)[4]
                
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
    
    def get_Zsensorsens(self):
        '''
        gets Zsensorsens for Z data.
        
        This is the sensitivity needed to convert the LSB data in nanometers for the Z-vs-T data chunk.
        '''        
        self.textfile.seek(0)
        
        zsensorsens_expr=re.compile(".*Sens. ZSensorSens")
        
        for line in self.textfile.readlines():
            if zsensorsens_expr.match(line):
                zsensorsens_raw_expression=line.split()
                #we must take only first occurrence, so we exit from the cycle immediately
                break
        
        return (float(zsensorsens_raw_expression[3]))*(10**(-9))
        
    def Z_data(self):
        '''
        returns converted ext and ret Z curves.
        They're on the second chunk (Z vs t).
        '''
        #Zmagnify_zt=self.get_Z_magnify(2)
        #Zscale_zt=self.get_Z_scale(2)
        Zlsb_zt=self.get_Z_scale_LSB(2)
        #rampsize_zt=self.get_ramp_size(2)
        #rampoffset_zt=self.get_ramp_offset(2)
        zsensorsens=self.get_Zsensorsens()
        
        '''
        The magic formula that converts the Z data is:
        
        meters = LSB * V_lsb_conversion_factor * ZSensorSens
        '''
        
        #z_curves=[item*Zlsb_zt*zsensorsens for item in self.data_chunks[1].pair['ext']],[item*Zlsb_zt*zsensorsens for item in self.data_chunks[1].pair['ret']]
        z_curves=[item*Zlsb_zt*zsensorsens for item in self.data_chunks[1].ext()],[item*Zlsb_zt*zsensorsens for item in self.data_chunks[1].ret()]                
        z_curves=[DataChunk(item) for item in z_curves]
        return z_curves
    
    def Z_extremes(self):
        '''
        returns the extremes of the Z values
        '''
        zcurves=self.Z_data()
        z_extremes={}
        z_extremes['ext']=zcurves[0][0],zcurves[0][-1]
        z_extremes['ret']=zcurves[1][0],zcurves[1][-1]
        
        return z_extremes
        
    def Z_step(self):
        '''
        returns the calculated step between the Z values
        '''
        zrange={}
        zpoints={}

        z_extremes=self.Z_extremes()
        
        zrange['ext']=abs(z_extremes['ext'][0]-z_extremes['ext'][1])
        zrange['ret']=abs(z_extremes['ret'][0]-z_extremes['ret'][1])
        
        #We must take 1 from the calculated zpoints, or when I use the arange function gives me a point more
        #with the step. That is, if I have 1000 points, and I use arange(start,stop,step), I have 1001 points...
        #For cleanness, solution should really be when using arange, but oh well...
        zpoints['ext']=len(self.Z_data()[0])-1
        zpoints['ret']=len(self.Z_data()[1])-1
        #this syntax must become coherent!!
        return (zrange['ext']/zpoints['ext']),(zrange['ret']/zpoints['ret'])
        
    def Z_domains(self):
        '''
        returns the Z domains on which to plot the force data.
        
        The Z domains are returned as a single long DataChunk() extended list. The extension and retraction part
        can be extracted using ext() and ret() methods.       
        '''   
        x1step=self.Z_step()[0]
        x2step=self.Z_step()[1]           
        
        try:
            xext=arange(self.Z_extremes()['ext'][0],self.Z_extremes()['ext'][1],-x1step)
            xret=arange(self.Z_extremes()['ret'][0],self.Z_extremes()['ret'][1],-x2step)
        except:
            xext=arange(0,1)
            xret=arange(0,1)
            print 'picoforce.py: Warning. xext, xret domains cannot be extracted.'
                
        if not (len(xext)==len(xret)):
            if self.debug:
                #print warning
                print "picoforce.py: Warning. Extension and retraction domains have different sizes."
                print "length extension: ", len(xext)
                print "length retraction: ", len(xret)
                print "You cannot trust the resulting curve."
                print "Until a solution is found, I substitute the ext domain with the ret domain. Sorry."
            xext=xret
        
        return DataChunk(xext.tolist()+xret.tolist())
        
    def Z_scan_size(self):
        return self.get_Zscan_V_size()*self.get_Z_scan_sensitivity()
        
    def Z_start(self):
        return self.get_Zscan_V_start()*self.get_Z_scan_sensitivity()
    
    def ramp_size(self,whichchunk):
        '''
        to be implemented if needed
        '''
        raise "Not implemented yet."
        
    
    def ramp_offset(self,whichchunk):
        '''
        to be implemented if needed
        '''
        raise "Not implemented yet."
    
    def detriggerize(self, forcext):
        '''
        Cuts away the trigger-induced s**t on the extension curve.
        DEPRECATED
        cutindex=2
        startvalue=forcext[0]
        
        for index in range(len(forcext)-1,2,-2):  
           if forcext[index]>startvalue:
                cutindex=index
           else:
                break

        return cutindex
        '''
        return 0
        
    def is_me(self):
        '''
        self-identification of file type magic
        '''
        curve_file=file(self.filepath)
        header=curve_file.read(30)
        curve_file.close()
        
        if header[2:17] == 'Force file list': #header of a picoforce file
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
        
                
        force=self.LSB_to_force()
        zdomain=self.Z_domains()
        
        samples=self._get_samples_line()
        #cutindex=0
        #cutindex=self.detriggerize(force.ext())
        
        main_plot=lhc.PlotObject()
        
        main_plot.vectors=[[zdomain.ext()[0:samples], force.ext()[0:samples]],[zdomain.ret()[0:samples], force.ret()[0:samples]]]
        main_plot.normalize_vectors()
        main_plot.units=['meters','newton']
        main_plot.destination=0
        main_plot.title=self.filepath
        
        
        return [main_plot]
