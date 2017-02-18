#!/usr/bin/env python

'''
Input check routines.

Copyright (C) 2008 Alberto Gomez-Casado (University of Twente).

This program is released under the GNU General Public License version 2.
'''

from types import *



def safeinput (message, valid=[]):
    '''
    friendlier frontend for alphainput and numinput
    valid should be a list of 0...n values
    '''

    #if possible values are not listed we just ask for any non-null input 
    if len(valid)==0:
        return alphainput(message, '',1,[])
    
    
    if len(valid)>0:
        #if valid values are string we use alphainput, if it is only one we take as default
        if type(valid[0]) is StringType:
            if len(valid)==1:
                return alphainput(message, valid[0], 0,[])
            else:
                return alphainput(message,'', 1,valid)
            
        #if valid values are numbers we use numinput
        if type(valid[0]) is IntType:
            if len(valid)==1:
                return numinput(message,valid[0],1,[])
            else:
                return numinput(message,'',1,valid)
    
    

def alphainput (message, default, repeat, valid):
    '''
    message: prompt for the user
    default: return value if user input was not correct (and repeat=0)
    repeat: keeps asking user till it gets a valid input
    valid: list of allowed answers, empty list for "anything"
    ''' 
    if default and not repeat:
        print 'Press [enter] for default: ('+str(default)+')'
    reply=raw_input(message)
    if len(valid)>0:
        if reply in valid: 
            return reply
        else:
            if repeat==1:
                while reply not in valid:
                    reply=raw_input('You should enter any of these: '+ str(valid) +'\n'+ message)
                return reply
            else:
                return default
    else:
        if len(reply)>0:
            return reply
        else:
            if not repeat:
                return default
            else:
                while len(reply)==0:
                    print 'Try again'
                    reply=raw_input(message)
                return reply

                    

def checkalphainput (test, default, valid):
    #useful when input was taken form command args
    if len(valid)>0:
        if test in valid: 
            return test
        else:
            return default
    else:
        #TODO: raise exception?
        if len(test)>0:
            return test
        else:
            return default


def numinput(message, default, repeat, limits):
    '''
    message: prompt for the user
    default: return value if user input was not correct (and repeat=0)
    repeat: keeps asking user till it gets a valid input
    limits: pair of values, input is checked to be between them, empty list for "any number"
    ''' 
    if default and not repeat:
        print 'Press [enter] for default: '+str(default)
        
    reply=raw_input(message)
    
    try:
        intreply=int(reply)
    except:
        intreply=None
              
    if len(limits)==2:
        high=int(limits.pop())
        low=int(limits.pop())
        if intreply>=low and intreply <= high:
            return intreply
        else:
            if repeat==1:
                while intreply<low or intreply>high :
                    reply=raw_input('You should enter values between: '+ str(low)+' and '+str(high) +'\n'+ message)
                    try:
                        intreply=int(reply)
                    except:
                        intreply=None
                return intreply
            else:
                return default
    else:
        if intreply!=None:
            return intreply
        else:
            if not repeat:
                return default
            else:
                while intreply==None:
                    print 'Try again'
                    reply=raw_input(message)
                    try:
                        intreply=int(reply)
                    except:
                        intreply=None
                return intreply

def checknuminput(test,default,limits):
    #useful when input was taken from command args
    if len(limits)==2:
        high=int(limits.pop())
        low=int(limits.pop())
        if test>=low and test <= high:
            return int(test)
        else:
            return default
    else:
        if len(test)>0:
            return int(test)
        else:
            return default

