#!/usr/bin/env python

'''
TUTORIAL PLUGIN FOR HOOKE

This plugin contains example commands to teach how to write an Hooke plugin, including description of main Hooke
internals.
(c)Massimo Sandal 2007
'''

import libhookecurve as lhc

import numpy as np

'''
SYNTAX OF DATA TYPE DECLARATION:
    type = type of object
    [ type ] = list containing objects of type
    {typekey:typearg} = dictionary with keys of type typekey and args of type typearg
    ( type ) = tuple containing objects of type
'''


class tutorialCommands:
    '''
    Here we define the class containing all the Hooke commands we want to define
    in the plugin.
    
    Notice the class name!!
    The syntax is filenameCommands. That is, if your plugin is pluggy.py, your class
    name is pluggyCommands.
    
    Otherwise, the class will be ignored by Hooke.
    ''' 
    
    def _plug_init(self):
        '''
        This is the plugin initialization.
        When Hooke starts and the plugin is loaded, this function is executed.
        If there is something you need to do when Hooke starts, code it in this function.
        '''
        print 'I am the Tutorial plugin initialization!'
        
        #Here we initialize a local configuration variable; see plotmanip_absvalue() for explanation.
        self.config['tutorial_absvalue']=0 
        pass
    
    def do_nothing(self,args):
        '''
        This is a boring but working example of an actual Hooke command.
        A Hooke command is a function of the xxxxCommands class, which is ALWAYS defined
        this way:
        
        def do_nameofcommand(self,args)
        
        *do_            is needed to make Hooke understand this function is a command
        *nameofcommand  is how the command will be called in the Hooke command line.
        *self           is, well, self
        *args           is ALWAYS needed (otherwise Hooke will crash executing the command). We will see
                        later what args is.
        
        Note that if you now start Hooke with this plugin activated and you type in the Hooke command
        line "help nothing" you will see this very text as output. So the help of a command is a
        string comment below the function definition, like this one.
        
        Commands usually return None.
        '''
        print 'I am a Hooke command. I do nothing.'
        
    def do_printargs(self,args):
        '''
        This command prints the args you give to it.
        args is always a string, that contains everything you write after the command.
        So if you issue "mycommand blah blah 12345" args is "blah blah 12345".
        
        Again, args is needed in the definition even if your command does not use it.
        '''
        print 'You gave me those args: '+args
        
    def help_tutorial(self):
        '''
        This is a help function. 
        If you want a help function for something that is not a command, you can write a help
        function like this. Calling "help tutorial" will execute this function.
        '''
        print 'You called help_tutorial()'
        
    def do_environment(self,args):
        '''
        This plugin contains a panoramic of the Hooke command line environment variables,
        and prints their current value.
        '''
        
        '''self.current_list
        TYPE: [ libhookecurve.HookeCurve ], len=variable
        contains the actual playlist of Hooke curve objects.
        Each HookeCurve object represents a reference to a data file.
        We will see later in detail how do they work.
        '''
        print 'current_list length:',len(self.current_list)
        print 'current_list 0th:',self.current_list[0]
        
        '''self.pointer
        TYPE: int
        contains the index of
        the current curve in the playlist
        '''
        print 'pointer: ',self.pointer
        
        '''self.current
        TYPE: libhookecurve.HookeCurve
        contains the current curve displayed.
        We will see later how it works.
        '''
        print 'current:',self.current
        
        '''self.plots
        TYPE: [ libhookecurve.PlotObject ], len=1,2
        contains the current default plots.
        Each PlotObject contains all info needed to display 
        the plot: apart from the data vectors, the title, destination
        etc.
        Usually self.plots[0] is the default topmost plot, self.plots[1] is the
        accessory bottom plot.
        '''
        print 'plots:',self.plots
        
        '''self.config
        TYPE: { string:anything }
        contains the current Hooke configuration variables, in form of a dictionary.
        '''
        print 'config:',self.config
        
        '''self.plotmanip
        TYPE: [ function ]
        Contains the ordered plot manipulation functions.
        These functions are called to modify the default plot by default before it is plotted.
        self.plots contains the plot passed through the plot manipulators.
        We will see it better later.
        *YOU SHOULD NEVER MODIFY THAT*
        '''
        print 'plotmanip: ',self.plotmanip
        
        '''self.drivers
        TYPE: [ class ]
        Contains the plot reading drivers.
        *YOU SHOULD NEVER MODIFY THAT*
        '''
        print 'drivers: ',self.drivers
        
        '''self.frame
        TYPE: wx.Frame
        Contains the wx Frame of the GUI.
        ***NEVER, EVER TOUCH THAT.***
        '''
        print 'frame: ',self.frame
        
        '''self.list_of_events
        TYPE: { string:wx.Event }
        Contains the wx.Events to communicate with the GUI.
        Usually not necessary to use it, unless you want
        to create a GUI plugin.
        '''
        print 'list of events:',self.list_of_events
        
        '''self.events_from_gui
        TYPE: Queue.Queue
        Contains the Queue where data from the GUI is put.
        Usually not necessary to use it, unless you want
        to create a GUI plugin.
        '''
        print 'events from gui:',self.events_from_gui
        
        '''self.playlist_saved
        TYPE: Int (0/1) ; Boolean
        Flag that tells if the playlist has been saved or not.
        '''
        print 'playlist saved:',self.playlist_saved
        
        '''self.playlist_name
        TYPE: string
        Name of current playlist
        '''
        print 'playlist name:',self.playlist_name
        
        '''self.notes_saved
        TYPE: Int (0/1) ; Boolean
        Flag that tells if the playlist has been saved or not.
        '''
        print 'notes saved:',self.notes_saved
        

    def do_myfirstplot(self,args):
        '''
        In this function, we see how to create a PlotObject and send it to the screen.
        ***Read the code of PlotObject in libhookecurve.py before!***.
        '''
        
        #We generate some interesting data to plot for this example.
        xdata1=np.arange(-5,5,0.1)
        xdata2=np.arange(-5,5,0.1)
        ydata1=[item**2 for item in xdata1]
        ydata2=[item**3 for item in xdata2]
        
        #Create the object.
        #The PlotObject class lives in the libhookecurve library.
        myplot=lhc.PlotObject()
        '''
        The *data* of the plot live in the .vectors list. 
        
        plot.vectors is a multidimensional array:
        plot.vectors[0]=set1
        plot.vectors[1]=set2
        plot.vectors[2]=sett3
        etc.
        
        2 curves in a x,y plot are:
        [[[x1],[y1]],[[x2],[y2]]]
        for example:
            x1          y1              x2         y2
        [[[1,2,3,4],[10,20,30,40]],[[3,6,9,12],[30,60,90,120]]]
        x1 = self.vectors[0][0]
        y1 = self.vectors[0][1]
        x2 = self.vectors[1][0]
        y2 = self.vectors[1][1]
        '''
        #Pour 0-th dataset into myplot: 
        myplot.add_set(xdata1,ydata1)
        
        #Pour 1-st dataset into myplot: 
        myplot.add_set(xdata2,ydata2)
        
        #Add units to x and y axes
        #units=[string, string]
        myplot.units=['x axis','y axis']
        
        #Where do we want the plot? 0=top, 1=bottom
        myplot.destination=1
        
        '''Send it to the GUI.
        Note that you *have* to encapsulate it into a list, so you
        have to send [myplot], not simply myplot.
        
        You can also send more two plots at once
        self.send_plot([plot1,plot2])
        '''
        self._send_plot([myplot])
        

    def do_myfirstscatter(self,args):
        '''
        How to draw a scatter plot.
        '''
        #We generate some interesting data to plot for this example.
        xdata1=np.arange(-5,5,1)
        xdata2=np.arange(-5,5,1)
        ydata1=[item**2 for item in xdata1]
        ydata2=[item**3 for item in xdata2]
        
        myplot=lhc.PlotObject()
        myplot.add_set(xdata1,ydata1)
        myplot.add_set(xdata2,ydata2)
        
        
        #Add units to x and y axes
        myplot.units=['x axis','y axis']
        
        #Where do we want the plot? 0=top, 1=bottom
        myplot.destination=1
        
        '''None=standard line plot
        'scatter'=scatter plot
        By default, the styles attribute is an empty list. If you
        want to define a scatter plot, you must define all other
        plots as None or 'scatter', depending on what you want.
        
        Here we define the second set to be plotted as scatter,
        and the first to be plotted as line.
        
        Here we define also the colors to be the default Matplotlib colors
        '''
        myplot.styles=[None,'scatter']
        myplot.colors=[None,None]
        self._send_plot([myplot])
        

    def do_clickaround(self,args):
        '''
        Here we click two points on the curve and take some parameters from the points
        we have clicked.
        '''
        
        '''
        points = self._measure_N_points(N=Int, whatset=Int)
        *N = number of points to measure(1...n)
        *whatset = data set to measure (0,1...n)
        *points = a list of ClickedPoint objects, one for each point requested
        '''
        points=self._measure_N_points(N=2,whatset=1)
        print 'You clicked the following points.'
        
        '''
        These are the absolute coordinates of the
        point clicked. 
        [float, float] = x,y
        '''
        print 'Absolute coordinates:'
        print points[0].absolute_coords
        print points[1].absolute_coords
        print
        
        '''
        These are the coordinates of the points
        clicked, remapped on the graph.
        Hooke looks at the graph point which X
        coordinate is next to the X coordinate of
        the point measured, and uses that point
        as the actual clicked point.
        [float, float] = x,y
        '''
        print 'Coordinates on the graph:'
        print points[0].graph_coords
        print points[1].graph_coords
        print
        
        '''
        These are the indexes of the clicked points
        on the dataset vector.
        '''
        print 'Index of points on the graph:'
        print points[0].index
        print points[1].index
        
        
    def help_thedifferentplots(self):
        '''
        The *three* different default plots you should be familiar with
        in Hooke.
        
        Each plot contains of course the respective data in their
        vectors attribute, so here you learn also which data access for
        each situation.
        '''
        print '''
        1. THE RAW, CURRENT PLOTS
        
        self.current
        ---
        Contains the current libhookecurve.HookeCurve container object.
        A HookeCurve object defines only two default attributes:
        
        * self.current.path = string
        The path of the current displayed curve
        
        * self.current.curve = libhookecurve.Driver
        The curve object. This is not only generated by the driver,
        this IS a driver instance in itself.
        This means that in self.current.curve you can access the
        specific driver APIs, if you know them.
        
        And defines only one method:
        * self.current.identify()
        Fills in the self.current.curve object.
        See in the cycling tutorial.
        
        *****
        The REAL curve data actually lives in:
        ---
        * self.current.curve.default_plots() = [ libhooke.PlotObject ]
        Contains the raw PlotObject-s, as "spitted out" by the driver, without any
        intervention.
        This is as close to the raw data as Hooke gets.
        
        One or two plots can be spit out; they are always enclosed in a list.
        *****
        
        Methods of self.current.curve are:
        ---
        
        * self.current.curve.is_me()
        (Used by identify() only.)
        
        * self.current.curve.close_all()
        Closes all driver open files; see the cycling tutorial.
        '''
        
        print '''
        2. THE PROCESSED, DEFAULT PLOT
        
        The plot that is spitted out by the driver is *not* the usual default plot
        that is displayed by calling "plot" at the Hooke prompt.
        
        This is because the raw, driver-generated plot is usually *processed* by so called
        *plot processing* functions. We will see in the tutorial how to define
        them. 
        
        For example, in force spectroscopy force curves, raw data are automatically corrected
        for deflection. Other data can be, say, filtered by default.
                
        The default plots are accessible in 
        self.plots = [ libhooke.PlotObject ]
        
        self.plots[0] is usually the topmost plot
        self.plots[1] is usually the bottom plot (if present)
        '''
        
        print '''
        3. THE PLOT DISPLAYED RIGHT NOW.
        
        Sometimes the plots you are displaying *right now* is different from the previous
        two. You may have a fit trace, you may have issued some command that spits out
        a custom plot and you want to rework that, whatever. 
        
        You can obtain in any moment the plot currently displayed by Hooke by issuing
        
        PlotObject = self._get_displayed_plot(dest)
        * dest = Int (0/1)
        dest=0 : top plot
        dest=1 : bottom plot
        '''
    
    
    def do_cycling(self,args):
        '''
        Here we cycle through our playlist and print some info on the curves we find.
        Cycling through the playlist needs a bit of care to avoid memory leaks and dangling
        open files...
        
        Look at the source code for more information.
        '''
        
        def things_when_cycling(item):
            '''
            We encapsulate here everything has to open the actual curve file.
            By doing it all here, we avoid to do acrobacies when deleting objects etc.
            in the main loop: we do the dirty stuff here.
            '''
            
            '''
            identify()
        
            This method looks for the correct driver in self.drivers to use;
            and puts the curve content in the .curve attribute.
            Basically, until identify() is called, the HookeCurve object
            is just an empty shell. When identify() is called (usually by
            the Hooke plot routine), the HookeCurve object is "filled" with
            the actual curve.
            '''
          
            item.identify(self.drivers)
            
            '''
            After the identify(), item.curve contains the curve, and item.curve.default_plots() behaves exactly like
            self.current.curve.default_plots() -but for the given item.
            '''
            itplot=item.curve.default_plots()
            
            print 'length of X1 vector:',len(itplot[0].vectors[0][0]) #just to show something
            
            '''
            The following three lines are a magic spell you HAVE to do
            before closing the function.
            (Otherwise you will be plagued by unpredicatable, system-dependent bugs.)
            '''
            item.curve.close_all() #Avoid open files dangling
            del item.curve #Avoid memory leaks
            del item #Just be paranoid. Don't ask.
            
            return
        
        
        c=0
        for item in self.current_list:
            print 'Looking at curve ',c,'of',len(self.current_list)
            things_when_cycling(item)
            c+=1
        
        return
        
            
        
    def plotmanip_absvalue(self, plot, current, customvalue=None):
        '''
        This function defines a PLOT MANIPULATOR.
        A plot manipulator is a function that takes a plot in input, does something to the plot
        and returns the modified plot in output.
        The function, once plugged, gets automatically called everytime self.plots is updated
        
        For example, in force spectroscopy force curves, raw data are automatically corrected
        for deflection. Other data can be, say, filtered by default.
        
        To create and activate a plot manipulator you have to:
            * Write a function (like this) which name starts with "plotmanip_" (just like commands
              start with "do_")
            * The function must support four arguments:
              self : (as usual)
              plot : a plot object
              current : (usually not used, deprecated)
              customvalue=None : a variable containing custom value(s) you need for your plot manipulators.
            * The function must return a plot object.
            * Add an entry in hooke.conf: if your function is "plotmanip_something" you will have
              to add <something/> in the plotmanips section: example
            
            <plotmanips>
                <detriggerize/>
                <correct/>
                <median/>
                <something/>        
            </plotmanips>
            
            Important: Plot manipulators are *in pipe*: each plot manipulator output becomes the input of the next one.
            The order in hooke.conf *is the order* in which plot manipulators are connected, so in the example above
            we have:
            self.current.curve.default_plots() --> detriggerize --> correct --> median --> something --> self.plots
        '''
        
        '''
        Here we see what is in a configuration variable to enable/disable the plot manipulator as user will using
        the Hooke "set" command.
        Typing "set tutorial_absvalue 0" disables the plot manipulator; typing "set tutorial_absvalue 1" will enable it.
        '''
        if not self.config['tutorial_absvalue']:
            return plot
        
        #We do something to the plot, for demonstration's sake
        #If we needed variables, we would have used customvalue.
        plot.vectors[0][1]=[abs(i) for i in plot.vectors[0][1]]
        plot.vectors[1][1]=[abs(i) for i in plot.vectors[1][1]]
        
        #Return the plot object.
        return plot
            
        
#TODO IN TUTORIAL:
#how to add lines to an existing plot!!
#peaks
#configuration files
#gui plugins
