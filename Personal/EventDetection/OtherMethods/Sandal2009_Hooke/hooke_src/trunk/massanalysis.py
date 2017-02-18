#!/usr/bin/env python

'''
massanalysis.py

Global analysis of force curves with various parameters

Requires:
libpeakspot.py
flatfilts.py
'''


import libpeakspot as lps
import libhookecurve as lhc
import libhooke as lh
import numpy as np

import csv

class massanalysisCommands:

    def _plug_init(self):
        self.mass_variables={}        
        self.interesting_variables=['curve','firstpeak_distance','lastpeak_distance','Npeaks','median_distance','mean_distance']
        self._clean_data()
        
    def _clean_data(self):
        for variable in self.interesting_variables:
            self.mass_variables[variable]=[]

    def peak_position_from_contact(self, item, locations):
        '''
        calculates X distance of a peak from the contact point
        '''
        item.identify(self.drivers)        
        
        real_positions=[]
        cut_index=self.find_contact_point()
        
        #we assume the first is the plot with the force curve
        plot=item.curve.default_plots()[0]
        xret=plot.vectors[1][0]
        
        start_x=xret[cut_index]
        
        real_positions=[abs((xret[index])-(start_x)) for index in locations] 
        #close all open files
        item.curve.close_all()
        #needed to avoid *big* memory leaks!
        del item.curve
        del item
        return real_positions

    def do_maplist(self,args):
        '''
        MAPLIST
        (flatfilts.py)
        ----
        pass
        '''
        self._clean_data() #if we recall it, clean previous data!
        min_deviation=self.convfilt_config['mindeviation']
          
        
        c=0
        for item in self.current_list:
            try:
                peak_location,peak_size=self.exec_has_peaks(item, min_deviation)
                real_positions=self.peak_position_from_contact(item, peak_location)
                
                self.mass_variables['Npeaks'].append(len(peak_location))
                
                if len(peak_location) > 1:
                    self.mass_variables['firstpeak_distance'].append(min(real_positions))
                    self.mass_variables['lastpeak_distance'].append(max(real_positions))
                    
                    distancepeaks=[]
                    for index in range(len(real_positions)-1):
                        distancepeaks.append(real_positions[index+1]-real_positions[index])
                else:
                    self.mass_variables['firstpeak_distance'].append(0)
                    self.mass_variables['lastpeak_distance'].append(0)    
                    
                if len(peak_location) > 2:
                    self.mass_variables['median_distance'].append(np.median(distancepeaks))
                    self.mass_variables['mean_distance'].append(np.mean(distancepeaks))    
                else:
                    self.mass_variables['median_distance'].append(0)
                    self.mass_variables['mean_distance'].append(0)   
                
                print 'curve',c
            except SyntaxError:
                print 'curve',c,'not mapped'
                pass
            
            c+=1

    def do_plotmap(self,args):
        '''
        '''
        args=args.split()
        if len(args)>1:
            x=self.mass_variables[args[0]]
            y=self.mass_variables[args[1]]
        else:
            print 'Give me two arguments between those:'
            print self.interesting_variables
            return
            
        scattermap=lhc.PlotObject()
        scattermap.vectors=[[]]
        scattermap.vectors[0].append(x)
        scattermap.vectors[0].append(y)
        
        scattermap.units=[args[0],args[1]]
        scattermap.styles=['scatter']
        scattermap.destination=1
        
        self._send_plot([scattermap])
        
    def do_savemaps(self,args):
        '''
        args=filename
        '''
        
        '''
        def csv_write_cols(data, f):
            
            #from Bruno Desthuillers on comp.lang.python
            
            writer = csv.writer(f)
            keys = data.keys()
            writer.writerow(dict(zip(keys,keys)))
            for row in zip(*data.values()):
                writer.writerow(dict(zip(keys, row))) 
        '''
        
        f=open(args,'wb')
        lh.csv_write_dictionary(f,self.mass_variables)
        f.close()
        
        