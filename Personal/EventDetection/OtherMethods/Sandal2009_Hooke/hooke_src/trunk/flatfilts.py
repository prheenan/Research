#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
FLATFILTS

Force spectroscopy curves filtering of flat curves
Licensed under the GNU LGPL version 2

Other plugin dependencies:
procplots.py (plot processing plugin)
'''
#from libhooke import WX_GOOD
#import wxversion
#wxversion.select(WX_GOOD)

import xml.dom.minidom

#import wx
import scipy
import numpy
from numpy import diff
import matplotlib.pyplot as plt
#import pickle

import libpeakspot as lps
import libhookecurve as lhc


class flatfiltsCommands:
    
    def _plug_init(self):
        #configurate convfilt variables
        convfilt_configurator=ConvfiltConfig()
        
        #different OSes have different path conventions
        if self.config['hookedir'][0]=='/':
            slash='/' #a Unix or Unix-like system
        else:
            slash='\\' #it's a drive letter, we assume it's Windows
        
        self.convfilt_config=convfilt_configurator.load_config(self.config['hookedir']+slash+'convfilt.conf')
    
    def do_flatfilt(self,args):
        '''
        FLATFILT
        (flatfilts.py)
        Filters out flat (featureless) curves of the current playlist,
        creating a playlist containing only the curves with potential
        features.
        ------------
        Syntax:
        flatfilt [min_npks min_deviation]

        min_npks = minmum number of points over the deviation
        (default=4)

        min_deviation = minimum signal/noise ratio
        (default=9)

        If called without arguments, it uses default values, that
        should work most of the times.
        '''
        median_filter=7
        min_npks=4
        min_deviation=9
        
        args=args.split(' ')
        if len(args) == 2:
            min_npks=int(args[0])
            min_deviation=int(args[1])
        else:
            pass
        
        print 'Processing playlist...'
        notflat_list=[]
        
        c=0
        for item in self.current_list:
            c+=1
                                   
            try:
                notflat=self.has_features(item, median_filter, min_npks, min_deviation)
                print 'Curve',item.path, 'is',c,'of',len(self.current_list),': features are ',notflat
            except:
                notflat=False
                print 'Curve',item.path, 'is',c,'of',len(self.current_list),': cannot be filtered. Probably unable to retrieve force data from corrupt file.'
            
            if notflat:
                item.features=notflat
                item.curve=None #empty the item object, to further avoid memory leak
                notflat_list.append(item)
        
        if len(notflat_list)==0:
            print 'Found nothing interesting. Check your playlist, could be a bug or criteria could be too much stringent'
            return
        else:
            print 'Found ',len(notflat_list),' potentially interesting curves'
            print 'Regenerating playlist...'
            self.pointer=0
            self.current_list=notflat_list
            self.current=self.current_list[self.pointer]
            self.do_plot(0)
                 
    def has_features(self,item,median_filter,min_npks,min_deviation):
        '''
        decides if a curve is flat enough to be rejected from analysis: it sees if there
        are at least min_npks points that are higher than min_deviation times the absolute value
        of noise.
   
        Algorithm original idea by Francesco Musiani, with my tweaks and corrections.
        '''
        retvalue=False
        
        item.identify(self.drivers)        
        #we assume the first is the plot with the force curve
        #do the median to better resolve features from noise
        flat_plot=self.plotmanip_median(item.curve.default_plots()[0], item, customvalue=median_filter)
        flat_vects=flat_plot.vectors 
        item.curve.close_all()
        #needed to avoid *big* memory leaks!
        del item.curve
        del item
        
        #absolute value of derivate        
        yretdiff=diff(flat_vects[1][1])
        yretdiff=[abs(value) for value in yretdiff]
        #average of derivate values
        diffmean=numpy.mean(yretdiff)
        yretdiff.sort()
        yretdiff.reverse()
        c_pks=0
        for value in yretdiff:
            if value/diffmean > min_deviation:
                c_pks+=1
            else:
                break
                    
        if c_pks>=min_npks:
            retvalue = c_pks
        
        del flat_plot, flat_vects, yretdiff
        
        return retvalue

    ################################################################
    #-----CONVFILT-------------------------------------------------
    #-----Convolution-based peak recognition and filtering.
    #Requires the libpeakspot.py library
    
    def has_peaks(self, plot, abs_devs=None, maxpeak=True, window=10, nocontact=False):
        '''
        Finds peak position in a force curve.
        FIXME: should be moved in libpeakspot.py
        '''
        if abs_devs==None:
            abs_devs=self.convfilt_config['mindeviation']
        
        
        xret=plot.vectors[1][0]
        yret=plot.vectors[1][1]
        #Calculate convolution.
        convoluted=lps.conv_dx(yret, self.convfilt_config['convolution'])
        
        #cut everything before the contact point
        cut_index=self.find_contact_point(plot)
	#with the curves without a contact region we don't want any cut
	if nocontact==True:
	  cut_index=0

        #cut even more, before the blind window
        start_x=xret[cut_index]
        blind_index=0
        for value in xret[cut_index:]:
            if abs((value) - (start_x)) > self.convfilt_config['blindwindow']*(10**-9):
                break
            blind_index+=1
        cut_index+=blind_index
        #do the dirty convolution-peak finding stuff
        noise_level=lps.noise_absdev(convoluted[cut_index:], self.convfilt_config['positive'], self.convfilt_config['maxcut'], self.convfilt_config['stable'])               
        above=lps.abovenoise(convoluted,noise_level,cut_index,abs_devs)     
        peak_location,peak_size=lps.find_peaks(above,seedouble=self.convfilt_config['seedouble'])
        #take the minimum or the maximum of a peak
        for i in range(len(peak_location)):
            peak=peak_location[i]
            valpk=min(yret[peak-window:peak+window])  #maximum in force (near the unfolding point)
            index_pk=yret[peak-window:peak+window].index(valpk)+(peak-window)            

            if maxpeak==False:
               valpk=max(yret[peak:peak+window]) #minimum in force, near the baseline
               index_pk=yret[peak:peak+window].index(valpk)+(peak)

#  Let's explain that for the minimum.  Immaging that we know that there is a peak at position/region 100 and you have found its y-value,
#  Now you look in the array, from 100-10 to 100+10  (if the window is 10).
#  This "100-10 to 100+10" is substancially a new array with its index. In this array you have 20
#  elements, so the index of your y-value will be 10.
#  Now to find the index in the TOTAL array you have to add the "position" of the "region" (that in this case
#  correspond to 100) and also substract the window size ---> (+100-10)

            peak_location[i]=index_pk
            
        return peak_location,peak_size
    
    
    def exec_has_peaks(self,item,abs_devs):
        '''
        encapsulates has_peaks for the purpose of correctly treating the curve objects in the convfilt loop,
        to avoid memory leaks
        '''
        item.identify(self.drivers)        
        #we assume the first is the plot with the force curve
        plot=item.curve.default_plots()[0]
        
        if 'flatten' in self.config['plotmanips']:
                    #If flatten is present, use it for better recognition of peaks...
                    flatten=self._find_plotmanip('flatten') #extract flatten plot manipulator
                    plot=flatten(plot, item, customvalue=1)
        
        peak_location,peak_size=self.has_peaks(plot,abs_devs)
        #close all open files
        item.curve.close_all()
        #needed to avoid *big* memory leaks!
        del item.curve
        del item
        return peak_location, peak_size
        
    #------------------------
    #------commands----------
    #------------------------    
    def do_peaks(self,args):
        '''
        PEAKS
        (flatfilts.py)
        Test command for convolution filter / test.
        ----
        Syntax: peaks [deviations]
        absolute deviation = number of times the convolution signal is above the noise absolute deviation.
        Default is 5.
        '''
        if len(args)==0:
            args=self.convfilt_config['mindeviation']
        
        try:
            abs_devs=float(args)
        except:
            print 'Wrong argument, using config value'
            abs_devs=float(self.convfilt_config['mindeviation'])
                        
        defplots=self.current.curve.default_plots()[0] #we need the raw, uncorrected plots
        
        if 'flatten' in self.config['plotmanips']:
            flatten=self._find_plotmanip('flatten') #extract flatten plot manipulator
            defplots=flatten(defplots, self.current)
        else:
            print 'You have the flatten plot manipulator not loaded. Enabling it could give you better results.'
        
        peak_location,peak_size=self.has_peaks(defplots,abs_devs)
        print 'Found '+str(len(peak_location))+' peaks.'
        to_dump='peaks '+self.current.path+' '+str(len(peak_location))
        self.outlet.push(to_dump)
        #print peak_location
        
        #if no peaks, we have nothing to plot. exit.
        if len(peak_location)==0:
            return
        
        #otherwise, we plot the peak locations.
        xplotted_ret=self.plots[0].vectors[1][0]
        yplotted_ret=self.plots[0].vectors[1][1]
        xgood=[xplotted_ret[index] for index in peak_location]
        ygood=[yplotted_ret[index] for index in peak_location]
        
        recplot=self._get_displayed_plot()
        recplot.vectors.append([xgood,ygood])
        if recplot.styles==[]:
            recplot.styles=[None,None,'scatter']
            recplot.colors=[None,None,None]
        else:
            recplot.styles+=['scatter']
            recplot.colors+=[None]
        
        self._send_plot([recplot])
        
    def do_convfilt(self,args):
        '''
        CONVFILT
        (flatfilts.py)
        Filters out flat (featureless) curves of the current playlist,
        creating a playlist containing only the curves with potential
        features.
        ------------
        Syntax:
        convfilt [min_npks min_deviation]

        min_npks = minmum number of peaks
        (to set the default, see convfilt.conf file; CONVCONF and SETCONF commands)

        min_deviation = minimum signal/noise ratio *in the convolution*
        (to set the default, see convfilt.conf file; CONVCONF and SETCONF commands)

        If called without arguments, it uses default values.
        '''
        
        min_npks=self.convfilt_config['minpeaks']
        min_deviation=self.convfilt_config['mindeviation']
        
        args=args.split(' ')
        if len(args) == 2:
            min_npks=int(args[0])
            min_deviation=int(args[1])
        else:
            pass
        
        print 'Processing playlist...'
        print '(Please wait)'
        notflat_list=[]
        
        c=0
        for item in self.current_list:
            c+=1
                                   
            try:    
                peak_location,peak_size=self.exec_has_peaks(item,min_deviation)
                if len(peak_location)>=min_npks:
                    isok='+'
                else:
                    isok=''
                print 'Curve',item.path, 'is',c,'of',len(self.current_list),': found '+str(len(peak_location))+' peaks.'+isok
            except:
                peak_location,peak_size=[],[]
                print 'Curve',item.path, 'is',c,'of',len(self.current_list),': cannot be filtered. Probably unable to retrieve force data from corrupt file.'
            
            if len(peak_location)>=min_npks:
                item.peak_location=peak_location
                item.peak_size=peak_size
                item.curve=None #empty the item object, to further avoid memory leak
                notflat_list.append(item)

            for i in range(1000):
		 k=0

        #Warn that no flattening had been done.
        if not ('flatten' in self.config['plotmanips']):
            print 'Flatten manipulator was not found. Processing was done without flattening.'
            print 'Try to enable it in your configuration file for better results.'
        
        if len(notflat_list)==0:
            print 'Found nothing interesting. Check your playlist, could be a bug or criteria could be too much stringent'
            return
        else:
            print 'Found ',len(notflat_list),' potentially interesting curves'
            print 'Regenerating playlist...'
            self.pointer=0
            self.current_list=notflat_list
            self.current=self.current_list[self.pointer]
            self.do_plot(0)
        
        
    def do_setconv(self,args):
        '''
        SETCONV
        (flatfilts.py)
        Sets the convfilt configuration variables
        ------
        Syntax: setconv variable value
        '''
        args=args.split()
        #FIXME: a general "set dictionary" function has to be built
        if len(args)==0:
            print self.convfilt_config
        else:
            if not (args[0] in self.convfilt_config.keys()):
                print 'This is not an internal convfilt variable!'
                print 'Run "setconv" without arguments to see a list of defined variables.'
                return
            
            if len(args)==1:
                print self.convfilt_config[args[0]]
            elif len(args)>1:
                try:
                    self.convfilt_config[args[0]]=eval(args[1])
                except NameError: #we have a string argument
                    self.convfilt_config[args[0]]=args[1]


#########################
#HANDLING OF CONFIGURATION FILE
class ConvfiltConfig:
    '''
    Handling of convfilt configuration file
    
    Mostly based on the simple-yet-useful examples of the Python Library Reference
    about xml.dom.minidom
    
    FIXME: starting to look a mess, should require refactoring
    '''
    
    def __init__(self):
        self.config={}
        
                
    def load_config(self, filename):
        myconfig=file(filename)                    
        #the following 3 lines are needed to strip newlines. otherwise, since newlines
        #are XML elements too, the parser would read them (and re-save them, multiplying
        #newlines...)
        #yes, I'm an XML n00b
        the_file=myconfig.read()
        the_file_lines=the_file.split('\n')
        the_file=''.join(the_file_lines)
                       
        self.config_tree=xml.dom.minidom.parseString(the_file)  
        
        def getText(nodelist):
            #take the text from a nodelist
            #from Python Library Reference 13.7.2
            rc = ''
            for node in nodelist:
                if node.nodeType == node.TEXT_NODE:
                    rc += node.data
            return rc
        
        def handleConfig(config):
            noiseabsdev_elements=config.getElementsByTagName("noise_absdev")
            convfilt_elements=config.getElementsByTagName("convfilt")
            handleAbsdev(noiseabsdev_elements)
            handleConvfilt(convfilt_elements)
                        
        def handleAbsdev(noiseabsdev_elements):
            for element in noiseabsdev_elements:
                for attribute in element.attributes.keys():
                    self.config[attribute]=element.getAttribute(attribute)
                    
        def handleConvfilt(convfilt_elements):
            for element in convfilt_elements:
                for attribute in element.attributes.keys():
                    self.config[attribute]=element.getAttribute(attribute)
            
        handleConfig(self.config_tree)
        #making items in the dictionary machine-readable
        for item in self.config.keys():
            try:
                self.config[item]=eval(self.config[item])
            except NameError: #if it's an unreadable string, keep it as a string
                pass
            
        return self.config
