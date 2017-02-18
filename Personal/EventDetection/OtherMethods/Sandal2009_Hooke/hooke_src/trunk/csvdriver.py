#!/usr/bin/env python

'''
csvdriver.py

Simple driver to read general comma-separated values in Hooke

Columns are read this way:
    
X1 , Y1 , X2 , Y2 , X3 , Y3 ...

If the number of columns is odd, the last column is ignored.

(c)Massimo Sandal, 2008
'''

import libhookecurve as lhc
import libhooke as lh
import csv

class csvdriverDriver(lhc.Driver):
    
        def __init__(self, filename):
        
            self.filedata = open(filename,'r')
            self.data = list(self.filedata) 
            self.filedata.close()
        
            self.filetype = 'generic'
            self.experiment = ''
            
            self.filename=filename
        
        def is_me(self):
            myfile=file(self.filename)
            headerline=myfile.readlines()[0]
            myfile.close()
            
            #using a custom header makes things much easier...
            #(looking for raw CSV data is at strong risk of confusion)
            if headerline[:-1]=='Hooke data':
                return True
            else:
                return False
        
        def close_all(self):
            self.filedata.close()
        
        def default_plots(self):
            rrows=csv.reader(self.data)
            rows=list(rrows) #transform the csv.reader iterator in a normal list
            columns=lh.transposed2(rows[1:])
            
            main_plot=lhc.PlotObject()
            main_plot.vectors=[]
            
            for index in range(0,len(columns),2):
                main_plot.vectors.append([])
                temp_x=columns[index]
                temp_y=columns[index+1]
                
                #convert to float (the csv gives strings)
                temp_x=[float(item) for item in temp_x]
                temp_y=[float(item) for item in temp_y]
                
                main_plot.vectors[-1].append(temp_x)
                main_plot.vectors[-1].append(temp_y)
                
            main_plot.units=['x','y']
            main_plot.title=self.filename
            main_plot.destination=0
            
            return [main_plot]
            
    