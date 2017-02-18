#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
hooke_cli.py

Command line module of Hooke.

Copyright (C) 2006 Massimo Sandal (University of Bologna, Italy).

This program is released under the GNU General Public License version 2.
'''


from libhooke import * #FIXME
import libhookecurve as lhc

import libinput as linp
import liboutlet as lout

from libhooke import WX_GOOD
from libhooke import HOOKE_VERSION

import wxversion
wxversion.select(WX_GOOD)
import wx

from wx.lib.newevent import NewEvent
from matplotlib.numerix import * #FIXME

import xml.dom.minidom
import sys, os, os.path, glob, shutil
import Queue
import cmd
import time

global __version__
global __codename__
global __releasedate__
__version__ = HOOKE_VERSION[0]
__codename__ = HOOKE_VERSION[1]
__releasedate__ = HOOKE_VERSION[2]

from matplotlib import __version__ as mpl_version
from wx import __version__ as wx_version
from wxmpl import __version__ as wxmpl_version
from scipy import __version__ as scipy_version
from numpy import __version__ as numpy_version
from sys import version as python_version
import platform


class HookeCli(cmd.Cmd):
    
    def __init__(self,frame,list_of_events,events_from_gui,config,drivers):
        cmd.Cmd.__init__(self)
                       
        self.prompt = 'hooke: '
        
        
        self.current_list=[] #the playlist we're using
        
        self.current=None    #the current curve under analysis. 
        self.plots=None
        '''
        The actual hierarchy of the "current curve" is a bit complex:
        
        self.current = the lhc.HookeCurve container object of the current curve
        self.current.curve = the current "real" curve object as defined in the filetype driver class
        self.current.curve.default_plots() = the default plots of the filetype driver.
        
        The plot objects obtained by mean of self.current.curve.default_plots() 
        then undergoes modifications by the plotmanip
        modifier functions. The modified plot is saved in self.plots and used if needed by other functions.       
        '''
        
        
        self.pointer=0       #a pointer to navigate the current list
                        
        #Things that come from outside
        self.frame=frame                        #the wx frame we refer to
        self.list_of_events=list_of_events      #a list of wx events we use to interact with the GUI
        self.events_from_gui=events_from_gui    #the Queue object we use to have messages from the GUI
        self.config=config                      #the configuration dictionary
        self.drivers=drivers                    #the file format drivers
        
        #get plot manipulation functions
        plotmanip_functions=[]
        for object_name in dir(self):
                if object_name[0:9]=='plotmanip':
                    plotmanip_functions.append(getattr(self,object_name))
        #put plotmanips in order
        self.plotmanip=[None for item in self.config['plotmanips']]
        for item in plotmanip_functions:
            namefunction=item.__name__[10:]
            if namefunction in self.config['plotmanips']:
                nameindex=self.config['plotmanips'].index(namefunction) #index of function in plotmanips config
                self.plotmanip[nameindex] = item
            else:
                pass
           
            
        self.playlist_saved=0 #Did we save the playlist?
        self.playlist_name='' #Name of playlist
        self.notes_saved=1 #Did we save the notes?
        self.notes_filename=None #Name of notes

        #create outlet
        self.outlet=lout.Outlet()
        
        #Data that must be saved in the playlist, related to the whole playlist (not individual curves)
        self.playlist_generics={} 
        
        #make sure we execute _plug_init() for every command line plugin we import
        for plugin_name in self.config['plugins']:
            try:
                plugin=__import__(plugin_name)
                try:
                    eval('plugin.'+plugin_name+'Commands._plug_init(self)')
                except AttributeError:
                    pass
            except ImportError:
                pass

        #load default list, if possible
        self.do_loadlist(self.config['defaultlist'])
        
#HELPER FUNCTIONS
#Everything sending an event should be here
    def _measure_N_points(self, N, whatset=1):
        '''
        general helper function for N-points measures
        '''
        wx.PostEvent(self.frame,self.list_of_events['measure_points'](num_of_points=N, set=whatset))
        while 1:
            try:
                points=self.frame.events_from_gui.get()
                break
            except Empty:
                pass
        return points
        
    def _get_displayed_plot(self,dest=0):
        '''
        returns the currently displayed plot.
        '''
        wx.PostEvent(self.frame, self.list_of_events['get_displayed_plot'](dest=dest))
        while 1:
            try:
                displayed_plot=self.events_from_gui.get()
            except Empty:
                pass
            if displayed_plot:
                break
        return displayed_plot
    
    def _send_plot(self,plots):
        '''
        sends a plot to the GUI
        '''
        wx.PostEvent(self.frame, self.list_of_events['plot_graph'](plots=plots))
        return
        
    def _find_plotmanip(self, name):
        '''
        returns a plot manipulator function from its name
        '''
        return self.plotmanip[self.config['plotmanips'].index(name)]
    
    def _clickize(self, xvector, yvector, index):
        '''
        returns a ClickedPoint() object from an index and vectors of x, y coordinates       
        '''
        point=ClickedPoint()
        point.index=index
        point.absolute_coords=xvector[index],yvector[index]
        point.find_graph_coords(xvector,yvector)
        return point
    
#HERE COMMANDS BEGIN
    
    def help_set(self):
        print '''
SET
Sets a local configuration variable
-------------
Syntax: set [variable] [value]
        '''
    def do_set(self,args):
        #FIXME: some variables in self.config should be hidden or intelligently configurated...
        args=args.split()
        if len(args)==0:
            print 'You must specify a variable and a value'
            print 'Available variables:'
            print self.config.keys()
            return
        if args[0] not in self.config.keys():
            print 'This is not an internal Hooke variable!'
            return
        if len(args)==1:
            #FIXME:we should reload the config file and reset the config value
            print self.config[args[0]]
            return
        key=args[0]
        try: #try to have a numeric value
            value=float(args[1])
        except ValueError: #if it cannot be converted to float, it's None, or a string...
            value=args[1]
            if value.lower()=='none':
                value=None
            else:
                value=args[1]
                
        self.config[key]=value
        self.do_plot(0)
        
#PLAYLIST MANAGEMENT AND NAVIGATION
#------------------------------------
    
    def help_loadlist(self):
        print '''
LOADLIST
Loads a file playlist
-----------
Syntax: loadlist [playlist file]
        '''
    def do_loadlist(self, args):
        #checking for args: if nothing is given as input, we warn and exit.
        while len(args)==0:
            args=linp.safeinput('File to load?')
        
        arglist=args.split()
        play_to_load=arglist[0]
        
        #We assume a Hooke playlist has the extension .hkp
        if play_to_load[-4:] != '.hkp':
            play_to_load+='.hkp'
        
        try:            
            playxml=PlaylistXML()
            self.current_list, self.playlist_generics=playxml.load(play_to_load)
            self.current_playxml=playxml
        except IOError:
            print 'File not found.'
            return
        
        print 'Loaded %s curves' %len(self.current_list)
        
        if 'pointer' in self.playlist_generics.keys():
            self.pointer=int(self.playlist_generics['pointer'])
        else:
            #if no pointer is found, set the current curve as the first curve of the loaded playlist
            self.pointer=0
        print 'Starting at curve ',self.pointer
            
        self.current=self.current_list[self.pointer]
        
        #resets saved/notes saved state
        self.playlist_saved=0
        self.playlist_name=''
        self.notes_saved=0        
    
        self.do_plot(0)
        
        
    def help_genlist(self):
        print '''
GENLIST
Generates a file playlist.
Note it doesn't *save* it: see savelist for this.

If [input files] is a directory, it will use all files in the directory for playlist.
So:
genlist dir
genlist dir/
genlist dir/*.*

are all equivalent syntax.
------------
Syntax: genlist [input files]
        
'''
    def do_genlist(self,args):
        #args list is: input path, output name
        if len(args)==0:
            args=linp.safeinput('Input files?')
                    
        arglist=args.split()      
        list_path=arglist[0]
                  
        #if it's a directory, is like /directory/*.*
        #FIXME: probably a bit kludgy.
        if os.path.isdir(list_path): 
            if platform.system == 'Windows':
                SLASH="\\"
            else:
                SLASH="/"
            if list_path[-1] == SLASH:
                list_path=list_path+'*'
            else:    
                list_path=list_path+SLASH+'*'
        
        #expanding correctly the input list with the glob module :)        
        list_files=glob.glob(list_path)
        list_files.sort()

        self.current_list=[]
        for item in list_files:
            try:
                if os.path.isfile(item):
                    self.current_list.append(lhc.HookeCurve(os.path.abspath(item))) 
            except:
                pass
            
        self.pointer=0    
        if len(self.current_list)>0:
            self.current=self.current_list[self.pointer]
        else:
            print 'Empty list!'
            return
        
        #resets saved/notes saved state
        self.playlist_saved=0
        self.playlist_name=''
        self.notes_saved=0  
        
        self.do_plot(0)
       
        
    def do_savelist(self,args):
        '''
        SAVELIST
        Saves the current file playlist on disk.
        ------------
        Syntax: savelist [filename]
        '''
        while len(args)==0:
            args=linp.safeinput('Output file?',['savedlist.txt'])
    
        output_filename=args
        
        self.playlist_generics['pointer']=self.pointer
        
        #autocomplete filename if not specified
        if output_filename[-4:] != '.hkp':
            output_filename+='.hkp'
        
        playxml=PlaylistXML()
        playxml.export(self.current_list, self.playlist_generics)
        playxml.save(output_filename)                  
        
        #remembers we have saved playlist
        self.playlist_saved=1
        
    def help_addtolist(self):
        print '''
ADDTOLIST
Adds a file to the current playlist
--------------
Syntax: addtolist [filename]
'''
    def do_addtolist(self,args):
        #args list is: input path
        if len(args)==0:
            print 'You must give the input filename you want to add'
            self.help_addtolist()
            return
          
        filenames=glob.glob(args)
        
        for filename in filenames:
            self.current_list.append(lhc.HookeCurve(os.path.abspath(filename)))
        #we need to save playlist
        self.playlist_saved=0
    
    def help_printlist(self):
        print '''
PRINTLIST
Prints the list of curves in the current playlist
-------------
Syntax: printlist
'''
    def do_printlist(self,args):
        for item in self.current_list:
            print item.path
            
    
    def help_jump(self):
        print '''
JUMP
Jumps to a given curve.
------
Syntax: jump {$curve}

If the curve is not in the current playlist, it politely asks if we want to add it.
        '''  
    def do_jump(self,filename):
        '''
        jumps to the curve with the given filename.
        if the filename is not in the playlist, it asks if we must add it or not.
        '''
        
        if filename=='':
            filename=linp.safeinput('Jump to?')
            
        filepath=os.path.abspath(filename)
        print filepath
                
        c=0
        item_not_found=1
        while item_not_found:
            try:
                
                if self.current_list[c].path == filepath:
                    self.pointer=c
                    self.current=self.current_list[self.pointer]
                    item_not_found=0
                    self.do_plot(0)
                else:
                    c+=1  
            except IndexError:
                #We've found the end of the list.
                answer=linp.safeinput('Curve not found in playlist. Add it to list?',['y'])
                if answer.lower()[0]=='y':
                    try:
                        self.do_addtolist(filepath)
                    except:
                        print 'Curve file not found.'
                        return
                    self.current=self.current_list[-1]
                    self.pointer=(len(current_list)-1)
                    self.do_plot(0)
                    
                item_not_found=0
    
    
    def do_index(self,args):
        '''
        INDEX
        Prints the index of the current curve in the list
        -----
        Syntax: index
        '''
        print self.pointer+1, 'of', len(self.current_list) 
    
    
    def help_next(self):
        print '''
NEXT
Go the next curve in the playlist.
If we are at the last curve, we come back to the first.
-----
Syntax: next, n
        '''
    def do_next(self,args):
        try:
            self.current.curve.close_all()
        except:
            print 'No curve file loaded, currently!'
            print 'This should not happen, report to http://code.google.com/p/hooke'
            return
        
        if self.pointer == (len(self.current_list)-1):
            self.pointer=0
            print 'Playlist finished; back to first curve.'
        else:
            self.pointer+=1
        
        self.current=self.current_list[self.pointer]
        self.do_plot(0)
        
    
    def help_n(self):
        self.help_next()

    def do_n(self,args):
	try:
          self.do_next(args)
        except:
	  print "Error in the playlist, have you correctly generated it?"
        
    def help_previous(self,args):
        print '''
PREVIOUS
Go to the previous curve in the playlist.
If we are at the first curve, we jump to the last.
-------
Syntax: previous, p
    '''
    def do_previous(self,args):
        try:
            self.current.curve.close_all()
        except:
            print 'No curve file loaded, currently!'
            print 'This should not happen, report to http://code.google.com/p/hooke'
            return
        if self.pointer == 0:
            self.pointer=(len(self.current_list)-1)
            print 'Start of playlist; jump to last curve.' 
        else:
            self.pointer-=1
            
        self.current=self.current_list[self.pointer]
        self.do_plot(args)
        
            
    def help_p(self):
        self.help_previous()
    def do_p(self,args):
        self.do_previous(args)

        
#PLOT INTERACTION COMMANDS
#-------------------------------    
    def help_plot(self):
        print '''
PLOT
Plots the current force curve
-------
Syntax: plot
        '''
    def do_plot(self,args):
        
        self.current.identify(self.drivers)
        self.plots=self.current.curve.default_plots()
        try:
            self.plots=self.current.curve.default_plots()
        except Exception, e:
            print 'Unexpected error occurred in do_plot().'
            print e
            return
            
        #apply the plotmanip functions eventually present
        nplots=len(self.plots)
        c=0
        while c<nplots:
            for function in self.plotmanip: #FIXME: something strange happens about self.plotmanip[0]
                self.plots[c]=function(self.plots[c], self.current)
                
            self.plots[c].xaxes=self.config['xaxes'] #FIXME: in the future, xaxes and yaxes should be set per-plot
            self.plots[c].yaxes=self.config['yaxes']
                
            c+=1

        self._send_plot(self.plots)
        
    def _delta(self, set=1):
        '''
        calculates the difference between two clicked points
        '''
        print 'Click two points'
        points=self._measure_N_points(N=2, whatset=set)
        dx=abs(points[0].graph_coords[0]-points[1].graph_coords[0])
        dy=abs(points[0].graph_coords[1]-points[1].graph_coords[1])
        unitx=self.plots[points[0].dest].units[0]
        unity=self.plots[points[0].dest].units[1]
        return dx,unitx,dy,unity
        
    def do_delta(self,args):
        '''
        DELTA
        
        Measures the delta X and delta Y between two points.
        ----
        Syntax: delta
        '''
        dx,unitx,dy,unity=self._delta()
        print str(dx)+' '+unitx
        print str(dy)+' '+unity
    
    def _point(self, set=1):
        '''calculates the coordinates of a single clicked point'''

        print 'Click one point'
        point=self._measure_N_points(N=1, whatset=set)
        
        x=point[0].graph_coords[0]
        y=point[0].graph_coords[1]
        unitx=self.plots[point[0].dest].units[0]
        unity=self.plots[point[0].dest].units[1]
        return x,unitx,y,unity
        
    def do_point(self,args):
        '''
        POINT
        
        Returns the coordinates of a point on the graph.
        ----
        Syntax: point
        '''
        x,unitx,y,unity=self._point()
        print str(x)+' '+unitx
        print str(y)+' '+unity
        to_dump='point '+self.current.path+' '+str(x)+' '+unitx+', '+str(y)+' '+unity
        self.outlet.push(to_dump)    
   
        
    def do_close(self,args=None):
        '''
        CLOSE
        Closes one of the two plots. If no arguments are given, the bottom plot is closed.
        ------
        Syntax: close [top,bottom]
        '''
        if args=='top':
            to_close=0
        elif args=='bottom':
            to_close=1
        else:
            to_close=1
        
        close_plot=self.list_of_events['close_plot']
        wx.PostEvent(self.frame, close_plot(to_close=to_close))
        
    def do_show(self,args=None):
        '''
        SHOW
        Shows both plots.
        ''' 
        show_plots=self.list_of_events['show_plots']
        wx.PostEvent(self.frame, show_plots())
       
        
    
    #PLOT EXPORT AND MANIPULATION COMMANDS
    def help_export(self):
        print '''
EXPORT
Saves the current plot as an image file
---------------
Syntax: export [filename] {plot to export}

The supported formats are PNG and EPS; the file extension of the filename is automatically recognized
and correctly exported. Resolution is (for now) fixed at 150 dpi.

If you have a multiple plot, the optional plot to export argument tells Hooke which plot you want to export. If 0, the top plot is exported. If 1, the bottom plot is exported (Exporting both plots is still to implement)
        '''
    def do_export(self,args):
        #FIXME: the bottom plot doesn't have the title
        
        dest=0
        
        if len(args)==0:
            #FIXME: We have to go into the libinput stuff and fix it, for now here's a dummy replacement...
            #name=linp.safeinput('Filename?',[self.current.path+'.png'])
            name=raw_input('Filename? ')
        else:
            args=args.split()
            name=args[0]
            if len(args) > 1:
                dest=int(args[1]) 
                
        export_image=self.list_of_events['export_image']
        wx.PostEvent(self.frame, export_image(name=name, dest=dest))
        
        
    def help_txt(self):
        print '''
TXT
Saves the current curve as a text file
Columns are, in order:
X1 , Y1 , X2 , Y2 , X3 , Y3 ...

-------------
Syntax: txt [filename] {plot to export} or
	txt [filename] all
	all  : To save all the curves in different windows in a single file.
        '''
    def do_txt(self,args):
        
        def transposed2(lists, defval=0):
            '''
            transposes a list of lists, i.e. from [[a,b,c],[x,y,z]] to [[a,x],[b,y],[c,z]] without losing
            elements
            (by Zoran Isailovski on the Python Cookbook online)
            '''
            if not lists: return []
            return map(lambda *row: [elem or defval for elem in row], *lists)
        
        whichplot=0
        args=args.split()
        if len(args)==0:
            filename=linp.safeinput('Filename?',[self.current.path+'.txt'])
        else:
            filename=linp.checkalphainput(args[0],self.current.path+'.txt',[])
            try:
		if args[1]=="all":
		  whichplot="all"
                else:
                  whichplot=int(args[1])
            except:
                pass
        
	if whichplot!="all":
	    try:
		outofplot=self.plots[whichplot].vectors
	    except:
		print "Plot index out of range."
		return 0
	    columns=[]     
	    for dataset in self.plots[whichplot].vectors:
		for i in range(0,len(dataset)): 
		    columns.append([])
		    for value in dataset[i]:
			#columns[-1].append(str(value*(10**9)))                   
			columns[-1].append(str(value))
	    rows=transposed2(columns, 'nan')
	    rows=[' , '.join(item) for item in rows]
	    text='\n'.join(rows)
	    
	    txtfile=open(filename,'w+')
	    #Save units of measure in header
	    txtfile.write('X:'+self.plots[whichplot].units[0]+'\n')
	    txtfile.write('Y:'+self.plots[whichplot].units[1]+'\n')
	    txtfile.write(text)
	    txtfile.close()

        else:
	  columns=[]
          for wp in range(len(self.plots)):     
	    for dataset in self.plots[wp].vectors:
		for i in range(0,len(dataset)): 
		    columns.append([])
		    for value in dataset[i]:
			#columns[-1].append(str(value*(10**9)))                   
			columns[-1].append(str(value))
	    rows=transposed2(columns, 'nan')
	    rows=[' , '.join(item) for item in rows]
	    text='\n'.join(rows)

	    txtfile=open(filename,'w+')
	    #Save units of measure in header
            for i in range(len(self.plots)):
	      txtfile.write('X:'+self.plots[i].units[0]+'\n')
	      txtfile.write('Y:'+self.plots[i].units[1]+'\n')
	    txtfile.write(text)
	    txtfile.close()
        
    
    #LOGGING, REPORTING, NOTETAKING
    

    def do_note_old(self,args):
        '''
        NOTE_OLD
        **deprecated**: Use note instead. Will be removed in 0.9
        
        Writes or displays a note about the current curve.
        If [anything] is empty, it displays the note, otherwise it adds a note.
        The note is then saved in the playlist if you issue a savelist command
        ---------------
        Syntax: note_old [anything]        

        '''
        if args=='':
            print self.current_list[self.pointer].notes
        else:
            #bypass UnicodeDecodeError troubles
            try:
                args=args.decode('ascii')
            except:
                args=args.decode('ascii','ignore')
                if len(args)==0:
                    args='?'
                    
            self.current_list[self.pointer].notes=args
        self.notes_saved=0
            
            
    def do_note(self,args):
        '''
        NOTE
        
        Writes or displays a note about the current curve.
        If [anything] is empty, it displays the note, otherwise it adds a note.
        The note is then saved in the playlist if you issue a savelist command.
        ---------------
        Syntax: note_old [anything]        

        '''
        if args=='':
            print self.current_list[self.pointer].notes
        else:
            if self.notes_filename == None:
		if not os.path.exists(os.path.realpath('output')):
		    os.mkdir('output')
                self.notes_filename=raw_input('Notebook filename? ')
		self.notes_filename=os.path.join(os.path.realpath('output'),self.notes_filename)
                title_line='Notes taken at '+time.asctime()+'\n'
                f=open(self.notes_filename,'a')
                f.write(title_line)
                f.close()
                
            #bypass UnicodeDecodeError troubles    
            try:
               args=args.decode('ascii')
            except:
               args=args.decode('ascii','ignore')
               if len(args)==0:
                   args='?'
            self.current_list[self.pointer].notes=args
            
            f=open(self.notes_filename,'a+')
            note_string=(self.current.path+'  |  '+self.current.notes+'\n')
            f.write(note_string)
            f.close()
                           
    def help_notelog(self):
        print '''
NOTELOG
Writes a log of the notes taken during the session for the current
playlist
--------------        
Syntax notelog [filename]
'''        
    def do_notelog(self,args):
        
        if len(args)==0:
            args=linp.safeinput('Notelog filename?',['notelog.txt'])
            
        note_lines='Notes taken at '+time.asctime()+'\n'
        for item in self.current_list:
            if len(item.notes)>0:
                #FIXME: log should be justified
                #FIXME: file path should be truncated...
                note_string=(item.path+'  |  '+item.notes+'\n')
                note_lines+=note_string
                
        try:
            f=open(args,'a+')
            f.write(note_lines)
            f.close()
        except IOError, (ErrorNumber, ErrorMessage):
            print 'Error: notes cannot be saved. Catched exception:'
            print ErrorMessage
        
        self.notes_saved=1

    def help_copylog(self):
        print '''
COPYLOG
Moves the annotated curves to another directory
-----------
Syntax copylog [directory]
        '''
    def do_copylog(self,args):
        
        if len(args)==0:
            args=linp.safeinput('Destination directory?')  #TODO default
        
        mydir=os.path.abspath(args)
        if not os.path.isdir(mydir):
            print 'Destination is not a directory.'
            return
        
        for item in self.current_list:
            if len(item.notes)>0:
                try:
                    shutil.copy(item.path, mydir)
                except (OSError, IOError):
                    print 'Cannot copy file. '+item.path+' Perhaps you gave me a wrong directory?'

#OUTLET management
#-----------------
    def do_outlet_show(self,args):
        '''OUTLET_SHOW
        ---------
        Shows current content of outlet with index for reference
        '''
        self.outlet.printbuf()

    def do_outlet_undo(self, args):
        '''OUTLET_UNDO
        ---------
        Eliminates last entry in outlet
        '''
        print 'Erasing last entry'
        self.outlet.pop()

    def do_outlet_delete(self, args):
        '''OUTLET_DELETE
        Eliminates a particular entry from outlet
        Syntax: outlet_delete n
        '''
        if len(args)==0:
            print 'Index needed!, use outlet_show to know it'
        else:
            self.outlet.delete(args)

#OS INTERACTION COMMANDS
#-----------------    
    def help_dir(self):
        print '''
DIR, LS
Lists the files in the directory
---------
Syntax: dir [path]
          ls  [path]
        '''
    def do_dir(self,args):
        
        if len(args)==0:
            args='*'
        print glob.glob(args)
        
    def help_ls(self):
        self.help_dir(self)
    def do_ls(self,args):
        self.do_dir(args)
        
    def help_pwd(self):
        print '''
PWD
Gives the current working directory.
------------
Syntax: pwd
        '''
    def do_pwd(self,args):
        print os.getcwd()         
    
    def help_cd(self):
        print '''
CD
Changes the current working directory
-----
Syntax: cd
        '''
    def do_cd(self,args):
        mypath=os.path.abspath(args)
        try:
            os.chdir(mypath)
        except OSError:
            print 'I cannot access that directory.'
    
    
    def help_system(self):
        print '''
SYSTEM
Executes a system command line and reports the output
-----
Syntax system [command line]
        '''
        pass
    def do_system(self,args):
        waste=os.system(args)           
            
    def do_debug(self,args):
        '''
        this is a dummy command where I put debugging things
        '''
        print self.config['plotmanips']
        pass
            
    def help_current(self):
        print '''
CURRENT
Prints the current curve path.
------
Syntax: current
        '''
    def do_current(self,args):
        print self.current.path
        
    def do_info(self,args):
        '''
        INFO
        ----
        Returns informations about the current curve.
        '''
        print 'Path: ',self.current.path
        print 'Experiment: ',self.current.curve.experiment
        print 'Filetype: ',self.current.curve.filetype
        for plot in self.current.curve.default_plots():
            for set in plot.vectors:
                lengths=[len(item) for item in set]
                print 'Data set size: ',lengths
        
    def do_version(self,args):
        '''
        VERSION
        ------
        Prints the current version and codename, plus library version. Useful for debugging.
        '''     
        print 'Hooke '+__version__+' ('+__codename__+')'
        print 'Released on: '+__releasedate__
        print '---'
        print 'Python version: '+python_version
        print 'WxPython version: '+wx_version
        print 'wxMPL version: '+wxmpl_version
        print 'Matplotlib version: '+mpl_version
        print 'SciPy version: '+scipy_version
        print 'NumPy version: '+numpy_version
        print '---'
        print 'Platform: '+str(platform.uname())
        print '---'
        print 'Loaded plugins:',self.config['loaded_plugins']
        
    def help_exit(self):
        print '''
EXIT, QUIT
Exits the program cleanly.
------
Syntax: exit
Syntax: quit
'''    
    def do_exit(self,args):
        we_exit='N'
        
        if (not self.playlist_saved) or (not self.notes_saved):
            we_exit=linp.safeinput('You did not save your playlist and/or notes. Exit?',['n'])
        else:
            we_exit=linp.safeinput('Exit?',['y'])
        
        if we_exit[0].upper()=='Y':
            wx.CallAfter(self.frame.Close)
            sys.exit(0)
        else:
            return
    
    def help_quit(self):
        self.help_exit()
    def do_quit(self,args):
        self.do_exit(args)





if __name__ == '__main__':
    mycli=HookeCli(0)
    mycli.cmdloop()
