#!/usr/bin/env python

'''
Basic Viewer and ascii saver example

Copyright (C) 2008 Alberto Gomez-Casado (University of Twente).

This program is released under the GNU General Public License version 2.
'''


import liboutlet as lout
import libinput as linput

class Viewer(object):
    source=[]
    data=[]
    dtype='all'
    action=[]  #alias to call the actual viewer function, makes it general
    

    def setdtype(self, dt):
        #determines wich type of data will be retrieved from outlet
        self.dtype=dt	

    def show(self):
        #TODO should print only data matching 'type'
        self.source.printbuf()

    def getdata(self):
        #retrieves data from outlet
        self.data=self.source.read_type(self.dtype)



class Ascii(Viewer):
#example viewer, it just retrieves data and writes it to a text file
#probably this should be in viewer.py?

	def __init__(self,outref):
		self.source=outref
		#tells the viewer which outlet has the data (so far only one in hooke)
		self.action=self.dump
		#this allows to call dump (or any other function, depending on the viewer) from the CLI using 'vwaction'

	def dump(self):
		#retrieves and saves data
		self.getdata()
		destination=linput.safeinput('Enter filename:',['results.txt'])
		destfile=open(destination,'w+')
		destfile.write('\n'.join(self.data))
		destfile.close()
	
