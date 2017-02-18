#!/usr/bin/env python

'''
generaltccd.py

General utilities for TCCD stuff
'''

class generaltccdCommands:
    
    def plotmanip_threshold(self, plot, current, customvalue=False):
        '''
        Cuts from the plot everything below the threshold.
        Set the threshold with "set tccd_threshold"
        '''
        
        if current.curve.experiment != 'smfluo':
            return plot
        
        if not self.config['tccd_threshold'] and (not customvalue):
            return plot
        
        if customvalue:
            thresh=customvalue
        else:
            thresh=self.config['tccd_threshold']
        
        for set in plot.vectors:
            newy=[]
            for value in set[1]:
                if abs(value) < thresh:
                    newy.append(0)
                else:
                    newy.append(value)
            
            set[1]=newy
                    
        return plot
                

    def plotmanip_coincident(self,plot,current, customvalue=False):
        '''
        Shows only coincident events
        '''
        if current.curve.experiment != 'smfluo':
            return plot
        
        if not self.config['tccd_coincident'] and (not customvalue):
            return plot
        
        newred=[]
        newblue=[]
        for index in range(len(plot.vectors[0][1])):
            if abs(plot.vectors[0][1][index])>self.config['tccd_threshold'] and abs(plot.vectors[1][1][index])>self.config['tccd_threshold']:
                newred.append(plot.vectors[0][1][index])
                newblue.append(plot.vectors[1][1][index])
            else:
                newred.append(0)
                newblue.append(0)
                
        plot.vectors[0][1]=newred
        plot.vectors[1][1]=newblue
     
        return plot