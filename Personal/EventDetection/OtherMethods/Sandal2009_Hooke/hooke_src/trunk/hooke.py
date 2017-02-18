#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
HOOKE - A force spectroscopy review & analysis tool

(C) 2008 Massimo Sandal

Copyright (C) 2008 Massimo Sandal (University of Bologna, Italy).

This program is released under the GNU General Public License version 2.
'''

from libhooke import HOOKE_VERSION
from libhooke import WX_GOOD

import os

import wxversion
wxversion.select(WX_GOOD)
import wx
import wxmpl
from wx.lib.newevent import NewEvent

import matplotlib.numerix as nx
import scipy as sp

from threading import *
import Queue

from hooke_cli import HookeCli
from libhooke import *
import libhookecurve as lhc

#import file versions, just to know with what we're working...
from hooke_cli import __version__ as hookecli_version

global __version__
global events_from_gui
global config
global CLI_PLUGINS
global GUI_PLUGINS
global LOADED_PLUGINS
global PLOTMANIP_PLUGINS
global FILE_DRIVERS

__version__=HOOKE_VERSION[0]
__release_name__=HOOKE_VERSION[1]

events_from_gui=Queue.Queue() #GUI ---> CLI COMMUNICATION

print 'Starting Hooke.'
#CONFIGURATION FILE PARSING
config_obj=HookeConfig()
config=config_obj.load_config('hooke.conf')

#IMPORTING PLUGINS

CLI_PLUGINS=[]
GUI_PLUGINS=[]
PLOTMANIP_PLUGINS=[]
LOADED_PLUGINS=[]

plugin_commands_namespaces=[]
plugin_gui_namespaces=[]
for plugin_name in config['plugins']:
    try:
        plugin=__import__(plugin_name)
        try:
            eval('CLI_PLUGINS.append(plugin.'+plugin_name+'Commands)') #take Command plugin classes
            plugin_commands_namespaces.append(dir(eval('plugin.'+plugin_name+'Commands')))
        except:
            pass
        try:
            eval('GUI_PLUGINS.append(plugin.'+plugin_name+'Gui)') #take Gui plugin classes
            plugin_gui_namespaces.append(dir(eval('plugin.'+plugin_name+'Gui')))
        except:
            pass
    except ImportError:
        print 'Cannot find plugin ',plugin_name
    else:
        LOADED_PLUGINS.append(plugin_name)
        print 'Imported plugin ',plugin_name

#eliminate names common to all namespaces
for i in range(len(plugin_commands_namespaces)):
    plugin_commands_namespaces[i]=[item for item in plugin_commands_namespaces[i] if (item != '__doc__' and item != '__module__' and item != '_plug_init')]
#check for conflicts in namespaces between plugins
#FIXME: only in commands now, because I don't have Gui plugins to check
#FIXME: how to check for plugin-defined variables (self.stuff) ??
plugin_commands_names=[]
whatplugin_defines=[]
plugin_gui_names=[]
for namespace,plugin_name in zip(plugin_commands_namespaces, config['plugins']):
    for item in namespace:
        if item in plugin_commands_names:
            i=plugin_commands_names.index(item) #we exploit the fact index gives the *first* occurrence of a name...
            print 'Error. Plugin ',plugin_name,' defines a function already defined by ',whatplugin_defines[i],'!'
            print 'This should not happen. Please disable one or both plugins and contact the plugin authors to solve the conflict.'
            print 'Hooke cannot continue.'
            exit()
        else:
            plugin_commands_names.append(item)
            whatplugin_defines.append(plugin_name)


config['loaded_plugins']=LOADED_PLUGINS #FIXME: kludge -this should be global but not in config!
#IMPORTING DRIVERS
#FIXME: code duplication
FILE_DRIVERS=[]
LOADED_DRIVERS=[]
for driver_name in config['drivers']:
    try:
        driver=__import__(driver_name)
        try:
            eval('FILE_DRIVERS.append(driver.'+driver_name+'Driver)')
        except:
            pass
    except ImportError:
        print 'Cannot find driver ',driver_name
    else:
        LOADED_DRIVERS.append(driver_name)
        print 'Imported driver ',driver_name
config['loaded_drivers']=LOADED_DRIVERS

#LIST OF CUSTOM WX EVENTS FOR CLI ---> GUI COMMUNICATION
#FIXME: do they need to be here?
list_of_events={}

plot_graph, EVT_PLOT = NewEvent()
list_of_events['plot_graph']=plot_graph

plot_contact, EVT_PLOT_CONTACT = NewEvent()
list_of_events['plot_contact']=plot_contact

measure_points, EVT_MEASURE_POINTS = NewEvent()
list_of_events['measure_points']=measure_points

export_image, EVT_EXPORT_IMAGE = NewEvent()
list_of_events['export_image']=export_image

close_plot, EVT_CLOSE_PLOT = NewEvent()
list_of_events['close_plot'] = close_plot

show_plots, EVT_SHOW_PLOTS = NewEvent()
list_of_events['show_plots'] = show_plots

get_displayed_plot, EVT_GET_DISPLAYED_PLOT = NewEvent()
list_of_events['get_displayed_plot'] = get_displayed_plot
#------------

class CliThread(Thread):

    def __init__(self,frame,list_of_events):
        Thread.__init__(self)

        #here we have to put temporary references to pass to the cli object.
        self.frame=frame
        self.list_of_events=list_of_events

        self.debug=0 #to be used in the future

    def run(self):
        print '\n\nThis is Hooke, version',__version__ , __release_name__
        print
        print '(c) Massimo Sandal & others, 2006-2008. Released under the GNU Lesser General Public License Version 3'
        print 'Hooke is Free software.'
        print '----'
        print ''

        def make_command_class(*bases):
            #FIXME: perhaps redundant
            return type(HookeCli)("HookeCliPlugged", bases + (HookeCli,), {})
        cli = make_command_class(*CLI_PLUGINS)(self.frame,self.list_of_events,events_from_gui,config,FILE_DRIVERS)
        cli.cmdloop()

'''
GUI CODE

FIXME: put it in a separate module in the future?
'''
class MainMenuBar(wx.MenuBar):
    '''
    Creates the menu bar
    '''
    def __init__(self):
        wx.MenuBar.__init__(self)
        '''the menu description. the key of the menu is XX&Menu, where XX is a number telling
        the order of the menus on the menubar.
        &Menu is the Menu text
        the corresponding argument is ('&Item', 'itemname'), where &Item is the item text and itemname
        the inner reference to use in the self.menu_items dictionary.

        See create_menus() to see how it works

        Note: the mechanism on page 124 of "wxPython in Action" is less awkward, maybe, but I want
        binding to be performed later. Perhaps I'm wrong :)
        ''' 

        self.menu_desc={'00&File':[('&Open playlist','openplaymenu'),('&Exit','exitmenu')], 
                        '01&Edit':[('&Export text...','exporttextmenu'),('&Export image...','exportimagemenu')],
                        '02&Help':[('&About Hooke','aboutmenu')]}
        self.create_menus()

    def create_menus(self):
        '''
        Smartish routine to create the menu from the self.menu_desc dictionary
        Hope it's a workable solution for the future.
        '''
        self.menus=[] #the menu objects to append to the menubar
        self.menu_items={} #the single menu items dictionary, to bind to events

        names=self.menu_desc.keys() #we gotta sort, because iterating keys goes in odd order
        names.sort()

        for name in names:
            self.menus.append(wx.Menu())
            for menu_item in self.menu_desc[name]:
                self.menu_items[menu_item[1]]=self.menus[-1].Append(-1, menu_item[0])

        for menu,name in zip(self.menus,names):
            self.Append(menu,name[2:])

class MainPanel(wx.Panel):
    def __init__(self,parent,id):  

        wx.Panel.__init__(self,parent,id)
        self.splitter = wx.SplitterWindow(self)

ID_FRAME=100        
class MainWindow(wx.Frame):
    '''we make a frame inheriting wx.Frame and setting up things on the init'''
    def __init__(self,parent,id,title):

        #-----------------------------
        #WX WIDGETS INITIALIZATION

        wx.Frame.__init__(self,parent,ID_FRAME,title,size=(800,600),style=wx.DEFAULT_FRAME_STYLE|wx.NO_FULL_REPAINT_ON_RESIZE)

        self.mainpanel=MainPanel(self,-1)
        self.cpanels=[]

        self.cpanels.append(wx.Panel(self.mainpanel.splitter,-1))
        self.cpanels.append(wx.Panel(self.mainpanel.splitter,-1))

        self.statusbar=wx.StatusBar(self,-1)
        self.SetStatusBar(self.statusbar)

        self.mainmenubar=MainMenuBar()
        self.SetMenuBar(self.mainmenubar)

        self.controls=[]
        self.figures=[]
        self.axes=[]

        #This is our matplotlib plot
        self.controls.append(wxmpl.PlotPanel(self.cpanels[0],-1))
        self.controls.append(wxmpl.PlotPanel(self.cpanels[1],-1))
        #These are our figure and axes, so to have easy references
        #Also, we initialize
        self.figures=[control.get_figure() for control in self.controls]
        self.axes=[figure.gca() for figure in self.figures]

	for i in range(len(self.axes)):
	  self.axes[i].xaxis.set_major_formatter(EngrFormatter())
	  self.axes[i].yaxis.set_major_formatter(EngrFormatter(2))


        self.cpanels[1].Hide()
        self.mainpanel.splitter.Initialize(self.cpanels[0])

        self.sizer_dance() #place/size the widgets

        self.controls[0].SetSize(self.cpanels[0].GetSize())
        self.controls[1].SetSize(self.cpanels[1].GetSize())

        #resize the frame to properly draw on Windows
        frameSize=self.GetSize()
        frameSize.DecBy(1, 1)
        self.SetSize(frameSize)
        '''
        #if you need the exact same size as before DecBy, uncomment this block
        frameSize.IncBy(1, 1)
        self.SetSize(frameSize)
        '''

        #-------------------------------------------
        #NON-WX WIDGETS INITIALIZATION

        #Flags.
        self.click_plot=0

        #FIXME: These could become a single flag with different (string?) values
        #self.on_measure_distance=False
        #self.on_measure_force=False

        self.plot_fit=False

        #Number of points to be clicked
        self.num_of_points = 2

        #Data.
        '''
            self.current_x_ext=[[],[]]
            self.current_y_ext=[[],[]]
            self.current_x_ret=[[],[]]
            self.current_y_ret=[[],[]]


            self.current_x_unit=[None,None]
            self.current_y_unit=[None,None]
            '''

        #Initialize xaxes, yaxes
        #FIXME: should come from config
        self.current_xaxes=0
        self.current_yaxes=0

        #Other


        self.index_buffer=[]

        self.clicked_points=[]

        self.measure_set=None

        self.events_from_gui = events_from_gui

        '''
            This dictionary keeps all the flags and the relative functon names that
            have to be called when a point is clicked.
            That is:
            - if point is clicked AND foo_flag=True
            - foo()

            Conversely, foo_flag is True if a corresponding event is launched by the CLI.

            self.ClickedPoints() takes care of handling this
            '''

        self.click_flags_functions={'measure_points':[False, 'MeasurePoints']}

        #Binding of custom events from CLI --> GUI functions!                       
        #FIXME: Should use the self.Bind() syntax
        EVT_PLOT(self, self.PlotCurve)
        EVT_PLOT_CONTACT(self, self.PlotContact)
        EVT_GET_DISPLAYED_PLOT(self, self.OnGetDisplayedPlot)
        EVT_MEASURE_POINTS(self, self.OnMeasurePoints)
        EVT_EXPORT_IMAGE(self,self.ExportImage)
        EVT_CLOSE_PLOT(self, self.OnClosePlot)
        EVT_SHOW_PLOTS(self, self.OnShowPlots)

        #This event and control decide what happens when I click on the plot 0.
        wxmpl.EVT_POINT(self, self.controls[0].GetId(), self.ClickPoint0)
        wxmpl.EVT_POINT(self, self.controls[1].GetId(), self.ClickPoint1)

        #RUN PLUGIN-SPECIFIC INITIALIZATION
        #make sure we execute _plug_init() for every command line plugin we import
        for plugin_name in config['plugins']:
            try:
                plugin=__import__(plugin_name)
                try:
                    eval('plugin.'+plugin_name+'Gui._plug_init(self)')
                    pass
                except AttributeError:
                    pass
            except ImportError:
                pass



    #WX-SPECIFIC FUNCTIONS
    def sizer_dance(self):
        '''
            adjust size and placement of wxpython widgets.
            '''
        self.splittersizer = wx.BoxSizer(wx.VERTICAL)
        self.splittersizer.Add(self.mainpanel.splitter, 1, wx.EXPAND)

        self.plot1sizer = wx.BoxSizer()
        self.plot1sizer.Add(self.controls[0], 1, wx.EXPAND)

        self.plot2sizer = wx.BoxSizer()
        self.plot2sizer.Add(self.controls[1], 1, wx.EXPAND)

        self.panelsizer=wx.BoxSizer()
        self.panelsizer.Add(self.mainpanel, -1, wx.EXPAND)

        self.cpanels[0].SetSizer(self.plot1sizer)
        self.cpanels[1].SetSizer(self.plot2sizer)

        self.mainpanel.SetSizer(self.splittersizer)
        self.SetSizer(self.panelsizer)

    def binding_dance(self):
        self.Bind(wx.EVT_MENU, self.OnOpenPlayMenu, self.menubar.menu_items['openplaymenu'])
        self.Bind(wx.EVT_MENU, self.OnExitMenu, self.menubar.menu_items['exitmenu'])
        self.Bind(wx.EVT_MENU, self.OnExportText, self.menubar.menu_items['exporttextmenu'])
        self.Bind(wx.EVT_MENU, self.OnExportImage, self.menubar.menu_items['exportimagemenu'])
        self.Bind(wx.EVT_MENU, self.OnAboutMenu, self.menubar.menu_items['aboutmenu'])

    # DOUBLE PLOT MANAGEMENT
    #----------------------
    def show_both(self):
        '''
            Shows both plots.
            '''
        self.mainpanel.splitter.SplitHorizontally(self.cpanels[0],self.cpanels[1])
        self.mainpanel.splitter.SetSashGravity(0.5)
        self.mainpanel.splitter.SetSashPosition(300) #FIXME: we should get it and restore it
        self.mainpanel.splitter.UpdateSize()

    def close_plot(self,plot):
        '''
            Closes one plot - only if it's open
            '''
        if not self.cpanels[plot].IsShown():
            return
        if plot != 0:
            self.current_plot_dest = 0
        else:
            self.current_plot_dest = 1
        self.cpanels[plot].Hide()
        self.mainpanel.splitter.Unsplit(self.cpanels[plot])
        self.mainpanel.splitter.UpdateSize()


    def OnClosePlot(self,event):
        self.close_plot(event.to_close)       

    def OnShowPlots(self,event):
        self.show_both()


    #FILE MENU FUNCTIONS
    #--------------------
    def OnOpenPlayMenu(self, event):
        pass 

    def OnExitMenu(self,event):
        pass

    def OnExportText(self,event):
        pass

    def OnExportImage(self,event):
        pass

    def OnAboutMenu(self,event):
        pass

    #PLOT INTERACTION    
    #----------------                        
    def PlotCurve(self,event):
        '''
            plots the current ext,ret curve.
            '''
        dest=0

        #FIXME: BAD kludge following. There should be a well made plot queue mechanism, with replacements etc.
        #---
        #If we have only one plot in the event, we already have one in self.plots and this is a secondary plot,
        #do not erase self.plots but append the new plot to it.
        if len(event.plots) == 1 and event.plots[0].destination != 0 and len(self.plots) == 1:
            self.plots.append(event.plots[0])
        #if we already have two plots and a new secondary plot comes, we substitute the previous
        if len(event.plots) == 1 and event.plots[0].destination != 0 and len(self.plots) > 1:
            self.plots[1] = event.plots[0]
        else:
            self.plots = event.plots

        #FIXME. Should be in PlotObject, somehow
        c=0
        for plot in self.plots:
            if self.plots[c].styles==[]:
                self.plots[c].styles=[None for item in plot.vectors] 
            if self.plots[c].colors==[]:
                self.plots[c].colors=[None for item in plot.vectors] 

        for plot in self.plots:
            '''
            MAIN LOOP FOR ALL PLOTS (now only 2 are allowed but...)
            '''
            if 'destination' in dir(plot):
                dest=plot.destination

            #if the requested panel is not shown, show it
            if not ( self.cpanels[dest].IsShown() ):
                self.show_both()

            self.axes[dest].hold(False)
            self.current_vectors=plot.vectors
            self.current_title=plot.title
            self.current_plot_dest=dest #let's try this way to take into account the destination plot...

            c=0

            if len(plot.colors)==0:
                plot.colors=[None] * len(plot.vectors)
            if len(plot.styles)==0:
                plot.styles=[None] * len(plot.vectors)     

            for vectors_to_plot in self.current_vectors: 
                if plot.styles[c]=='scatter':
                    if plot.colors[c]==None:
                        self.axes[dest].scatter(vectors_to_plot[0], vectors_to_plot[1])
                    else:
                        self.axes[dest].scatter(vectors_to_plot[0], vectors_to_plot[1],color=plot.colors[c])
                else:
                    if plot.colors[c]==None:
                        self.axes[dest].plot(vectors_to_plot[0], vectors_to_plot[1])
                    else:
                        self.axes[dest].plot(vectors_to_plot[0], vectors_to_plot[1], color=plot.colors[c])
                self.axes[dest].hold(True)
                c+=1

            '''
                for vectors_to_plot in self.current_vectors:
                    if len(vectors_to_plot)==2: #3d plots are to come...
                        if len(plot.styles) > 0 and plot.styles[c] == 'scatter':
                            self.axes[dest].scatter(vectors_to_plot[0],vectors_to_plot[1])
                        elif len(plot.styles) > 0 and plot.styles[c] == 'scatter_red':
                            self.axes[dest].scatter(vectors_to_plot[0],vectors_to_plot[1],color='red')
                        else:
                            self.axes[dest].plot(vectors_to_plot[0],vectors_to_plot[1])

                        self.axes[dest].hold(True)
                        c+=1
                    else:
                        pass
                '''               
            #FIXME: tackles only 2d plots
            self.axes[dest].set_xlabel(plot.units[0])
            self.axes[dest].set_ylabel(plot.units[1])

            #FIXME: set smaller fonts
            self.axes[dest].set_title(plot.title)

            if plot.xaxes: 
                #swap X axis
                xlim=self.axes[dest].get_xlim()
                self.axes[dest].set_xlim((xlim[1],xlim[0])) 
            if plot.yaxes:
                #swap Y axis
                ylim=self.axes[dest].get_ylim()        
                self.axes[dest].set_ylim((ylim[1],ylim[0])) 

	    for i in range(len(self.axes)):
	      self.axes[i].xaxis.set_major_formatter(EngrFormatter())
	      self.axes[i].yaxis.set_major_formatter(EngrFormatter(2))


            self.controls[dest].draw()


    def PlotContact(self,event):
        '''
            plots the contact point
            DEPRECATED!
            '''
        self.axes[0].hold(True)
        self.current_contact_index=event.contact_index

        #now we fake a clicked point 
        self.clicked_points.append(ClickedPoint())
        self.clicked_points[-1].absolute_coords=self.current_x_ret[dest][self.current_contact_index], self.current_y_ret[dest][self.current_contact_index]
        self.clicked_points[-1].is_marker=True    

        self._replot()
        self.clicked_points=[]

    def OnMeasurePoints(self,event):
        '''
            trigger flags to measure N points
            '''
        self.click_flags_functions['measure_points'][0]=True
        if 'num_of_points' in dir(event):
            self.num_of_points=event.num_of_points
        if 'set' in dir(event):    
            self.measure_set=event.set            

    def ClickPoint0(self,event):
        self.current_plot_dest=0
        self.ClickPoint(event)
    def ClickPoint1(self,event):
        self.current_plot_dest=1
        self.ClickPoint(event)

    def ClickPoint(self,event):
        '''
            this function decides what to do when we receive a left click on the axes.
            We trigger other functions:
            - the action chosen by the CLI sends an event
            - the event raises a flag : self.click_flags_functions['foo'][0]
            - the raised flag wants the function in self.click_flags_functions[1] to be called after a click
            '''
        for key, value in self.click_flags_functions.items():
            if value[0]:
                eval('self.'+value[1]+'(event)')



    def MeasurePoints(self,event,current_set=1):
        dest=self.current_plot_dest
        try:
            current_set=self.measure_set
        except AttributeError:
            pass

        #find the current plot matching the clicked destination
        plot=self._plot_of_dest()
        if len(plot.vectors)-1 < current_set: #what happens if current_set is 1 and we have only 1 vector?
            current_set=current_set-len(plot.vectors)

        xvector=plot.vectors[current_set][0]
        yvector=plot.vectors[current_set][1]

        self.clicked_points.append(ClickedPoint())            
        self.clicked_points[-1].absolute_coords=event.xdata, event.ydata
        self.clicked_points[-1].find_graph_coords(xvector,yvector)
        self.clicked_points[-1].is_marker=True    
        self.clicked_points[-1].is_line_edge=True
        self.clicked_points[-1].dest=dest                

        self._replot()

        if len(self.clicked_points)==self.num_of_points:
            self.events_from_gui.put(self.clicked_points)
            #restore to default state:
            self.clicked_points=[]
            self.click_flags_functions['measure_points'][0]=False    


    def OnGetDisplayedPlot(self,event):
        if 'dest' in dir(event):
            self.GetDisplayedPlot(event.dest)
        else:
            self.GetDisplayedPlot(self.current_plot_dest)

    def GetDisplayedPlot(self,dest):
        '''
            returns to the CLI the currently displayed plot for the given destination
            '''
        displayed_plot=self._plot_of_dest(dest)
        events_from_gui.put(displayed_plot)

    def ExportImage(self,event):
        '''
            exports an image as a file.
            Current supported file formats: png, eps
            (matplotlib docs say that jpeg should be supported too, but with .jpg it doesn't work for me!)
            '''
        #dest=self.current_plot_dest
        dest=event.dest
        filename=event.name
        self.figures[dest].savefig(filename)

    '''
        def _find_nearest_point(self, mypoint, dataset=1):

            #Given a clicked point on the plot, finds the nearest point in the dataset (in X) that
            #corresponds to the clicked point.

            dest=self.current_plot_dest

            xvector=plot.vectors[dataset][0]
            yvector=plot.vectors[dataset][1]

            #Ye Olde sorting algorithm...
            #FIXME: is there a better solution?
            index=0
            best_index=0
            best_diff=10^9 #hope we never go over this magic number :(
            for point in xvector:
                diff=abs(point-mypoint)
                if diff<best_diff:
                    best_index=index
                    best_diff=diff
                index+=1

            return best_index,xvector[best_index],yvector[best_index]
         '''   

    def _plot_of_dest(self,dest=None):
        '''
            returns the plot that has the current destination
            '''
        if dest==None:
            dest=self.current_plot_dest
        try:
          plot=None
          for aplot in self.plots:
              if aplot.destination == dest:
                  plot=aplot
          return plot
        except:
           print "No curve available"
           return None

    def _replot(self):
        '''
            this routine is needed for a fresh clean-and-replot of interface
            otherwise, refreshing works very badly :(

            thanks to Ken McIvor, wxmpl author!
            '''
        dest=self.current_plot_dest
        #we get current zoom limits
        xlim=self.axes[dest].get_xlim()
        ylim=self.axes[dest].get_ylim()           
        #clear axes
        self.axes[dest].cla()

        #Plot curve:         
        #find the current plot matching the clicked destination
        plot=self._plot_of_dest()
        #plot all superimposed plots 
        c=0 
        if len(plot.colors)==0:
            plot.colors=[None] * len(plot.vectors)
        if len(plot.styles)==0:
            plot.styles=[None] * len(plot.vectors)     
        for plotset in plot.vectors: 
            if plot.styles[c]=='scatter':
                if plot.colors[c]==None:
                    self.axes[dest].scatter(plotset[0], plotset[1])
                else:
                    self.axes[dest].scatter(plotset[0], plotset[1],color=plot.colors[c])
            else:
                if plot.colors[c]==None:
                    self.axes[dest].plot(plotset[0], plotset[1])
                else:
                    self.axes[dest].plot(plotset[0], plotset[1], color=plot.colors[c])
            '''    
                if len(plot.styles) > 0 and plot.styles[c]=='scatter':
                    self.axes[dest].scatter(plotset[0], plotset[1],color=plot.colors[c])
                elif len(plot.styles) > 0 and plot.styles[c] == 'scatter_red':
                    self.axes[dest].scatter(plotset[0],plotset[1],color='red')
                else:
                    self.axes[dest].plot(plotset[0], plotset[1])
                '''
            c+=1
        #plot points we have clicked
        for item in self.clicked_points:
            if item.is_marker:
                if item.graph_coords==(None,None): #if we have no graph coords, we display absolute coords
                    self.axes[dest].scatter([item.absolute_coords[0]],[item.absolute_coords[1]])
                else:
                    self.axes[dest].scatter([item.graph_coords[0]],[item.graph_coords[1]])               

        if self.plot_fit:
            print 'DEBUGGING WARNING: use of self.plot_fit is deprecated!'
            self.axes[dest].plot(self.plot_fit[0],self.plot_fit[1])

        self.axes[dest].hold(True)      
        #set old axes again
        self.axes[dest].set_xlim(xlim)
        self.axes[dest].set_ylim(ylim)
        #set title and names again...
        self.axes[dest].set_title(self.current_title)           
        self.axes[dest].set_xlabel(plot.units[0])
        self.axes[dest].set_ylabel(plot.units[1])
        #and redraw!
        self.controls[dest].draw()


class MySplashScreen(wx.SplashScreen):
    """
    Create a splash screen widget.
    That's just a fancy addition... every serious application has a splash screen!
    """
    def __init__(self, frame):
        # This is a recipe to a the screen.
        # Modify the following variables as necessary.
        #aBitmap = wx.Image(name = "wxPyWiki.jpg").ConvertToBitmap()
        aBitmap=wx.Image(name='hooke.jpg').ConvertToBitmap()
        splashStyle = wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT
        splashDuration = 2000 # milliseconds
        splashCallback = None
        # Call the constructor with the above arguments in exactly the
        # following order.
        wx.SplashScreen.__init__(self, aBitmap, splashStyle,
                                 splashDuration, None, -1)
        wx.EVT_CLOSE(self, self.OnExit)
        self.frame=frame
        wx.Yield()

    def OnExit(self, evt):
        self.Hide()

        self.frame.Show()
        # The program will freeze without this line.
        evt.Skip()  # Make sure the default handler runs too...


#------------------------------------------------------------------------------

def main():

    #save the directory where Hooke is located
    config['hookedir']=os.getcwd()

    #now change to the working directory.
    try:
        os.chdir(config['workdir'])
    except OSError:
        print "Warning: Invalid work directory."

    app=wx.PySimpleApp()

    def make_gui_class(*bases):
        return type(MainWindow)("MainWindowPlugged", bases + (MainWindow,), {})

    main_frame = make_gui_class(*GUI_PLUGINS)(None, -1, ('Hooke '+__version__))

    #FIXME. The frame.Show() is called by the splashscreen here! Ugly as hell.

    mysplash=MySplashScreen(main_frame)
    mysplash.Show()

    my_cmdline=CliThread(main_frame, list_of_events)
    my_cmdline.start()


    app.MainLoop()

main()
