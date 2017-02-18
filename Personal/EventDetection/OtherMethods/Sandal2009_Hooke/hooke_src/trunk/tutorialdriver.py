#!/usr/bin/env python

'''
tutorialdriver.py

TUTORIAL DRIVER FOR HOOKE

Example driver to teach how to write a driver for data types.
(c)Massimo Sandal 2008
'''

'''
Here we define a (fake) file format that is read by this driver. The file format is as following:

TUTORIAL_FILE
PLOT1
X1
n1
n2
...
nN
Y1
n1
n2
...
nN
X2
n1
n2
..
nN
Y2
n1
n2
..
nN
PLOT2
X1
...
Y1
...
X2
...
Y2
...
END
that is, two plots with two datasets each.
'''

import libhookecurve as lhc #We need to import this library to define some essential data types

class tutorialdriverDriver(lhc.Driver):
    '''
    This is a *generic* driver, not a specific force spectroscopy driver.
    See the written documentation to see what a force spectroscopy driver must be defined to take into account Hooke facilities.
    
    Our driver is a class with the name convention nameofthedriverDriver, where "nameofthedriver" is the filename.
    The driver must inherit from the parent class lhc.Driver, so the syntax is
    class nameofthedriverDriver(lhc.Driver)
    '''
    def __init__(self, filename):
        '''
        THIS METHOD MUST BE DEFINED.
        The __init__ method MUST call the filename, so that it can open the file.
        '''
        self.filename=filename #self.filename can always be useful, and should be defined
        self.filedata = open(filename,'r') #We open the file
        '''
        In this case, we have a data format that is just a list of ASCII values, so we can just divide that in rows, and generate a list
        with each item being a row.
        Of course if your data files are binary, or follow a different approach, do whatever you need. :) 
        '''
        self.data = list(self.filedata)
        self.filedata.close() #remember to close the file
        
        '''These are two strings that can be used by Hooke commands/plugins to understand what they are looking at. They have no other
        meaning. They have to be somehow defined however - commands often look for those variables.
        
        self.filetype should contain the name of the exact filetype defined by the driver (so that filetype-specific commands can know
                      if they're dealing with the correct filetype)
        self.experiment should contain instead the type of data involved (for example, various drivers can be used for force-clamp experiments,
                      but hooke commands could like to know if we're looking at force clamp data, regardless of their origin, and not other 
                      kinds of data)
        
        Of course, all other variables you like can be defined in the class.
        '''
        self.filetype = 'tutorial'
        self.experiment = 'generic'
    
    def is_me(self):
        '''
        THIS METHOD MUST BE DEFINED.
        RETURNS: Boolean (True or False)
        This method must be an heuristic that looks at the file content and decides if the file can be opened by the driver itself. 
        It returns True if the file opened can be interpreted by the current driver, False otherwise. 
        Defining this method allows Hooke to understand what kind of files we're looking at automatically.
        
        We have to re-open/re-close the file here.
        '''
        
        myfile=open(self.filename, 'r')
        headerline=myfile.readlines()[0] #we take the first line
        myfile.close()
            
        '''
        Here, our "magic fingerprint" is the TUTORIAL_FILE header. Of course, depending on the data file, you can have interesting
        headers, or patterns, etc. that you can use to guess the data format. What matters is successful recognizing, and returning
        a boolean (True/False).
        '''
        if headerline[:-1]=='TUTORIAL_FILE': #[:-1], otherwise the return character is included in the line
            return True
        else:
            return False
        
    def _generate_vectors(self):
        '''
        Here we parse the data and generate the raw vectors. This method has only to do with the peculiar file format here, so it's of
        no big interest (I just wrote it to present a functional driver). 
        
        Only thing to remember, it can be nice to call methods that are used only "internally" by the driver (or by plugins) with a
        "_" prefix, so to have a visual remark. But it's just an advice.
        '''
        vectors={'PLOT1':[[],[],[],[]] , 'PLOT2':[[],[],[],[]]}
        positions={'X1':0,'Y1':1,'X2':2,'Y2':3}
        whatplot=''
        pos=0
        for item in self.data:
            try:
                num=float(item)
                vectors[whatplot][pos].append(num)
            except ValueError:
                if item[:-1]=='PLOT1':
                    whatplot=item[:-1]
                elif item[:-1]=='PLOT2':
                    whatplot=item[:-1]
                elif item[0]=='X' or item[0]=='Y':
                    pos=positions[item[:-1]]         
                else:
                    pass
        
        return vectors

    def close_all(self):
        '''
        THIS METHOD MUST BE DEFINED.
        This method is a precaution method that is invoked when cycling to avoid eventually dangling open files.
        In this method, all file objects defined in the driver must be closed.
        '''
        self.filename.close()
    
    
    def default_plots(self):
        '''
        THIS METHOD MUST BE DEFINED.
        RETURNS: [ lhc.PlotObject ] or [ lhc.PlotObject, lhc.PlotObject]
        
        This is the method that returns the plots to Hooke.
        It must return a list with *one* or *two* PlotObjects.
        
        See the libhookecurve.py source code to see how PlotObjects are defined and work in detail.
        '''
        gen_vectors=self._generate_vectors()
        
        #Here we create the main plot PlotObject and initialize its vectors.
        main_plot=lhc.PlotObject()
        main_plot.vectors=[]
        #The same for the other plot.
        other_plot=lhc.PlotObject()
        other_plot.vectors=[]
        
        '''
        Now we fill the plot vectors with our data.
                                                           set 1                           set 2
        The "correct" shape of the vector is [ [[x1,x2,x3...],[y1,y2,y3...]] , [[x1,x2,x3...],[y1,y2,y3...]] ], so we have to put stuff in this way into it.
        
        The add_set() method takes care of this , just use plot.add_set(x,y).
        '''
        main_plot.add_set(gen_vectors['PLOT1'][0],gen_vectors['PLOT1'][1])
        main_plot.add_set(gen_vectors['PLOT1'][2],gen_vectors['PLOT1'][3])
        
        other_plot.add_set(gen_vectors['PLOT2'][0],gen_vectors['PLOT2'][1])
        other_plot.add_set(gen_vectors['PLOT2'][2],gen_vectors['PLOT2'][3])
                
        '''
        normalize_vectors() trims the vectors, so that if two x/y couples are of different lengths, the latest
        points are trimmed (otherwise we have a python error). Always a good idea to run it, to avoid crashes.
        '''
        main_plot.normalize_vectors()
        other_plot.normalize_vectors()
        
        '''
        Here we define:
        - units: [string, string], define the measure units of X and Y axes
        - destination: 0/1 , defines where to plot the plot (0=top, 1=bottom), default=0
        - title: string , the plot title.
        
        for each plot.
        Again, see libhookecurve.py comments for details.
        '''
        main_plot.units=['unit of x','unit of y']
        main_plot.destination=0
        main_plot.title=self.filename+' main'
        
        other_plot.units=['unit of x','unit of y']
        other_plot.destination=1
        other_plot.title=self.filename+' other'
        
        return [main_plot, other_plot]