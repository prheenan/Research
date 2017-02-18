#!/usr/bin/env python

'''
mfp3d.py

Driver for MFP-3D files.

Copyright 2010 by Dr. Rolf Schmidt (Concordia University, Canada)
This driver is based on the work of R. Naud and A. Seeholzer (see below)
to read Igor binary waves. Code used with permission.

Modified for usage with Hooke CLI by Alberto Gomez-Casado (University of Twente, The Netherlands)

This program is released under the GNU General Public License version 2.
'''

# DEFINITION:
# Reads Igor's (Wavemetric) binary wave format, .ibw, files.
#
# ALGORITHM:
# Parsing proper to version 2, 3, or version 5 (see Technical notes TN003.ifn:
# http://mirror.optus.net.au/pub/wavemetrics/IgorPro/Technical_Notes/) and data
# type 2 or 4 (non complex, single or double precision vector, real values).
#
# AUTHORS:
# Matlab version: R. Naud August 2008 (http://lcn.epfl.ch/~naud/Home.html)
# Python port: A. Seeholzer October 2008
#
# VERSION: 0.1
#
# COMMENTS:
# Only tested for version 2 Igor files for now, testing for 3 and 5 remains to be done.
# More header data could be passed back if wished. For significance of ignored bytes see
# the technical notes linked above.

import numpy
import os.path
import struct

import libhookecurve as lhc


__version__='0.0.0.20100310'


class DataChunk(list):
    #Dummy class to provide ext and ret methods to the data list.
    
    def ext(self):
        halflen=(len(self)/2)
        return self[0:halflen]
        
    def ret(self):
        halflen=(len(self)/2)
        return self[halflen:]

class mfp3dDriver(lhc.Driver):

    #Construction and other special methods
    
    def __init__(self,filename):
        '''
        constructor method
        '''
           
	self.textfile    =file(filename)
        self.binfile=file(filename,'rb')
	#unnecesary, but some other part of the program expects these to be open     

        self.forcechunk=0
        self.distancechunk=1
	#TODO eliminate the need to set chunk numbers
        
        self.filepath=filename
        self.debug=True
        
	self.data = []
        self.note = []
        self.retract_velocity = None
        self.spring_constant = None
        self.filename = filename

        self.filedata = open(filename,'rU')
        self.lines = list(self.filedata.readlines())
        self.filedata.close()

        self.filetype = 'mfp3d'
        self.experiment = 'smfs'
             
     
    def _get_data_chunk(self,whichchunk):

	data = None
        f = open(self.filename, 'rb')
        ####################### ORDERING
        # machine format for IEEE floating point with big-endian
        # byte ordering
        # MacIgor use the Motorola big-endian 'b'
        # WinIgor use Intel little-endian 'l'
        # If the first byte in the file is non-zero, then the file is a WinIgor
        firstbyte = struct.unpack('b', f.read(1))[0]
        if firstbyte == 0:
            format = '>'
        else:
            format = '<'
        #######################  CHECK VERSION
        f.seek(0)
        version = struct.unpack(format+'h', f.read(2))[0]
        #######################  READ DATA AND ACCOMPANYING INFO
        if version == 2 or version == 3:
            # pre header
            wfmSize = struct.unpack(format+'i', f.read(4))[0] # The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.
            noteSize = struct.unpack(format+'i', f.read(4))[0] # The size of the note text.
            if version==3:
                formulaSize = struct.unpack(format+'i', f.read(4))[0]
            pictSize = struct.unpack(format+'i', f.read(4))[0] # Reserved. Write zero. Ignore on read.
            checksum = struct.unpack(format+'H', f.read(2))[0] # Checksum over this header and the wave header.
            # wave header
            dtype = struct.unpack(format+'h', f.read(2))[0]
            if dtype == 2:
                dtype = numpy.float32(.0).dtype
            elif dtype == 4:
                dtype = numpy.double(.0).dtype
            else:
                assert False, "Wave is of type '%i', not supported" % dtype
            dtype = dtype.newbyteorder(format)

            ignore = f.read(4) # 1 uint32
            bname = self._flatten(struct.unpack(format+'20c', f.read(20)))
            ignore = f.read(4) # 2 int16
            ignore = f.read(4) # 1 uint32
            dUnits = self._flatten(struct.unpack(format+'4c', f.read(4)))
            xUnits = self._flatten(struct.unpack(format+'4c', f.read(4)))
            npnts = struct.unpack(format+'i', f.read(4))[0]
            amod = struct.unpack(format+'h', f.read(2))[0]
            dx = struct.unpack(format+'d', f.read(8))[0]
            x0 = struct.unpack(format+'d', f.read(8))[0]
            ignore = f.read(4) # 2 int16
            fsValid = struct.unpack(format+'h', f.read(2))[0]
            topFullScale = struct.unpack(format+'d', f.read(8))[0]
            botFullScale = struct.unpack(format+'d', f.read(8))[0]
            ignore = f.read(16) # 16 int8
            modDate = struct.unpack(format+'I', f.read(4))[0]
            ignore = f.read(4) # 1 uint32
            # Numpy algorithm works a lot faster than struct.unpack
            data = numpy.fromfile(f, dtype, npnts)

        elif version == 5:
            # pre header
            checksum = struct.unpack(format+'H', f.read(2))[0] # Checksum over this header and the wave header.
            wfmSize = struct.unpack(format+'i', f.read(4))[0] # The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.
            formulaSize = struct.unpack(format+'i', f.read(4))[0]
            noteSize = struct.unpack(format+'i', f.read(4))[0] # The size of the note text.
            dataEUnitsSize = struct.unpack(format+'i', f.read(4))[0]
            dimEUnitsSize = struct.unpack(format+'4i', f.read(16))
            dimLabelsSize = struct.unpack(format+'4i', f.read(16))
            sIndicesSize = struct.unpack(format+'i', f.read(4))[0]
            optionSize1 = struct.unpack(format+'i', f.read(4))[0]
            optionSize2 = struct.unpack(format+'i', f.read(4))[0]

            # header
            ignore = f.read(4)
            CreationDate =  struct.unpack(format+'I',f.read(4))[0]
            modData =  struct.unpack(format+'I',f.read(4))[0]
            npnts =  struct.unpack(format+'i',f.read(4))[0]
            # wave header
            dtype = struct.unpack(format+'h',f.read(2))[0]
            if dtype == 2:
                dtype = numpy.float32(.0).dtype
            elif dtype == 4:
                dtype = numpy.double(.0).dtype
            else:
                assert False, "Wave is of type '%i', not supported" % dtype
            dtype = dtype.newbyteorder(format)

            ignore = f.read(2) # 1 int16
            ignore = f.read(6) # 6 schar, SCHAR = SIGNED CHAR?         ignore = fread(fid,6,'schar'); #
            ignore = f.read(2) # 1 int16
            bname = self._flatten(struct.unpack(format+'32c',f.read(32)))
            ignore = f.read(4) # 1 int32
            ignore = f.read(4) # 1 int32
            ndims = struct.unpack(format+'4i',f.read(16)) # Number of of items in a dimension -- 0 means no data.
            sfA = struct.unpack(format+'4d',f.read(32))
            sfB = struct.unpack(format+'4d',f.read(32))
            dUnits = self._flatten(struct.unpack(format+'4c',f.read(4)))
            xUnits = self._flatten(struct.unpack(format+'16c',f.read(16)))
            fsValid = struct.unpack(format+'h',f.read(2))
            whpad3 = struct.unpack(format+'h',f.read(2))
            ignore = f.read(16) # 2 double
            ignore = f.read(40) # 10 int32
            ignore = f.read(64) # 16 int32
            ignore = f.read(6) # 3 int16
            ignore = f.read(2) # 2 char
            ignore = f.read(4) # 1 int32
            ignore = f.read(4) # 2 int16
            ignore = f.read(4) # 1 int32
            ignore = f.read(8) # 2 int32

            data = numpy.fromfile(f, dtype, npnts)
            note_str = f.read(noteSize)
            note_lines = note_str.split('\r')
            self.note = {}
            for line in note_lines:
            	if ':' in line:
            	    key, value = line.split(':', 1)
                    self.note[key] = value
            self.retract_velocity = float(self.note['Velocity'])
            self.spring_constant = float(self.note['SpringConstant'])
        else:
            assert False, "Fileversion is of type '%i', not supported" % dtype
            data = []

        f.close()
        if len(data) > 0:
            #we have 3 columns: deflection, LVDT, raw
            #TODO detect which is each one
	    count = npnts / 3
            lvdt = data[:count] 
            deflection = data[count:2 * count] 
            #every column contains data for extension and retraction
            #we assume the same number of points for each
            #we could possibly extract this info from the note
            count = npnts / 6

	    forcechunk=deflection*self.spring_constant
	    distancechunk=lvdt
	
	    if	whichchunk==self.forcechunk:
		return forcechunk
	    if whichchunk==self.distancechunk:
		return distancechunk
        else:
            return None                          
        
    def _force(self):
	#returns force vector
        Kspring=self.spring_constant
        return DataChunk([(meter*Kspring) for meter in self._deflection()])

    def _deflection(self):
        #for internal use (feeds _force)
        deflect=self.data_chunks[self.forcechunk]/self.spring_constant      	     
        return deflect

    def _flatten(self, tup):
        out = ''
        for ch in tup:
            out += ch
        return out            
    
    def _Z(self):   
        return DataChunk(self.data_chunks[self.distancechunk])
        
    def is_me(self):
        if len(self.lines) < 34:
            return False

        name, extension = os.path.splitext(self.filename)
        if extension == '.ibw':
            for line in self.lines:
                if line.startswith('ForceNote:'):
		    self.data_chunks=[self._get_data_chunk(num) for num in [0,1,2]]
                    return True
            else:
                return False
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
        main_plot=lhc.PlotObject()
        main_plot.vectors=[[zdomain.ext(), force.ext()],[zdomain.ret(), force.ret()]]
        main_plot.normalize_vectors()
        main_plot.units=['meters','newton']
        main_plot.destination=0
        main_plot.title=self.filepath
        
        
        return [main_plot]

    def deflection(self):
        #interface for correct plotmanip and others
        deflectionchunk=DataChunk(self._deflection())
        return deflectionchunk.ext(),deflectionchunk.ret()
