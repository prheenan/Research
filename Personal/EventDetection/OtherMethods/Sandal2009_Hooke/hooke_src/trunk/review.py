#!/usr/bin/env python

'''review.py
Alberto Gomez-Casado (c) 2010 University of Twente
'''

from libhooke import WX_GOOD
import wxversion
wxversion.select(WX_GOOD)
from wx import PostEvent
import numpy as np
import shutil
import libinput as linp
import copy
import os.path

import warnings
warnings.simplefilter('ignore',np.RankWarning)


class reviewCommands:

    def do_review(self,args):
        '''
        REVIEW
        (review.py)
        Presents curves (in current playlist) in groups of ten. User can indicate which curves will be selected to be saved in a separate directory for further analysis.
	By default curves are presented separated -30 nm in x and -100 pN in y. 
	Curve number one of each set is the one showing the approach.
        ------------
        Syntax:
        review [x spacing (nm)] [y spacing (pN]

        '''

	args=args.split()

	if len(args)==2:
		try:
			xgap=int(args[0])*1e-9 #scale to SI units 
			ygap=int(args[1])*1e-12
		except:
			print 'Spacings must be numeric! Using defaults'
			xgap=-30*1e-9
			ygap=-100*1e-12 
	else:
		xgap=-30*1e-9
		ygap=-100*1e-12 
	          
        print 'Processing playlist...'
        print '(Please wait)'
        keep_list=[]
        
	c=0
	print 'You can stop the review at any moment by entering \'q\' you can go back ten curves entering \'b\''
	print 'What curve no. you would like to start? (Enter for starting from the first)'
	skip=raw_input()
	
	if skip.isdigit()==False:
	    skip=0
	else:
	    skip=int(skip)
	    print 'Skipping '+str(skip)+ ' curves'
	    c=skip	

        while c < len(self.current_list):
             	
	
	    #take a group of ten curves and plot them with some spacing
	    			
	    curveset=self.current_list[c:c+10]
	
	    base=curveset[0]	
	    self.current=base
	    self.do_plot(0)	
	    multiplot=copy.deepcopy(self._get_displayed_plot(0))
	    self.current.curve.close_all()	

	    for i in range(1,10):
		if i >= len(curveset):
			print 'End of the list'
			print 'WARNING: maybe you want to finish!'
			break
		nextitem=curveset[i]
		if not nextitem.identify(self.drivers):
			continue		
		nextplot=self.plotmanip_correct(nextitem.curve.default_plots()[0],nextitem)
		nextvect=nextplot.vectors
		nextitem.curve.close_all()

		nextx=nextvect[1][0]
		nexty=nextvect[1][1]
		#center y around 0	
		ymedian=np.median(nexty)
		pos=0		
		for j in range(0,len(nextx)):
			nextx[j]=nextx[j]+i*xgap
			nexty[j]=nexty[j]+i*ygap-ymedian
		multiplot.add_set(nextx,nexty) 
        	multiplot.styles.append('lines')
        	multiplot.colors.append(None)
        	
	    self._send_plot([multiplot])		
	    			
	    
	    print 'Which ones you want to keep?'	
	    keep=raw_input()
	    if keep.isalpha():
		if keep=='b':
			print 'Going back ten curves'
			c-=10
			if c<0:
				print 'We are already at the start'
				c=0
			continue
		if keep=='q':
			break
	    else:
		for i in keep.split():
			if i.isdigit() and int(i)>0 and int(i)<11: #if it is not digit the int() call is never made, so no exception should be happening
				keep_item=curveset[int(i)-1].path
				if keep_item in keep_list:
					print 'This curve ('+keep_item+') was already selected, skipping'
				else:
					keep_list.append(keep_item)
			else:
				print 'You entered an invalid value: '+i
					
	    c+=10

	#FIXME I don't know why the print below gives errors sometimes
	try:
		print 'Kept '+str(len(keep_list))+' curves from '+str(min(c+i+1,len(self.current_list)))
	except:
		print 'Display error, never mind, we continue. Below the amount of kept and total curves:'
		print str(len(keep_list))
		print str(len(self.current_list))

	allok=0  #flag to keep from losing all the work in a slight mistake
	while allok==0:
		if len(keep_list) == 0:
			return
		save=linp.safeinput('Do you want to save the selected curves?',['y','n'])
		if save=='y':
			savedir=linp.safeinput('Destination directory?')
	        	savedir=os.path.abspath(savedir)
        		if not os.path.isdir(savedir):
            			print 'Destination is not a directory. Try again'
				continue
            	if save=='n':
			allok=1
			return
		
        	for item in keep_list:
            		try:
                    		shutil.copy(item, savedir)
				allok=1
                	except (OSError, IOError):
                    		print 'Cannot copy file. '+item+' Perhaps you gave me a wrong directory?'
				allok=0
				break

	return
