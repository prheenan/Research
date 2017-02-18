#!/usr/bin/env python

'''
Basic outlet object

Copyright (C) 2008 Alberto Gomez-Casado (University of Twente).

This program is released under the GNU General Public License version 2.
'''


import re


class Outlet(object):
    
    def __init__(self):
        self.buffer=[]
	#relations is still unused
        self.relations=[]
    
    def push(self, args):
	#queue new entry
        self.buffer.append(args)

    def pop(self):
	#delete last entry
	return self.buffer.pop();

    def printbuf(self):
	j=1;
        for i in self.buffer:
            print j, i
            j=j+1

    def delete(self, number):
	#delete entry matching given index
	if len(self.buffer)>int(number)-1 and int(number)>0:
	    self.buffer.pop(int(number)-1)		

    def empty(self):
        self.buffer=[]
        
    def read_last(self):
        return self.buffer[len(self.buffer)-1]
    
    def read_first(self):
        return self.buffer[0]

    def read_type(self,dtype):
	#returns entries matching a given type (force, distance, point...)        
	aux=[]
        index=0
	if dtype=='all':
		return self.buffer
	for i in self.buffer:
		if re.match(dtype+'*',i):
			aux.append(i)
        return aux
    
    
