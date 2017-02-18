# -*- coding: utf-8 -*-
class cutCommands:

    def _plug_init(self):
        self.cut_basecurrent=None
        self.cut_basepoints=None




    def do_cut(self,args):
        '''
CUT
        (cut.py)
        Cut the selected signal between two points.
	The first parameters is useful to select the window, with a single window wichplot is "0"  (zero).
        With the second parameter you have to select the signal (for FS for example
        you can select with "0" the approacing curve and 1 for the retracting
	curve. This depend also on how many set of data you have on the graph).
        With the second parameter you select the output name file for the selection.
	The data is arranged in two simple column without a header, the first column
	is the "x" data and the second the "y".
        -----------------
        Syntax: cut "whichplot" "whatset" "namefile"
        '''
        if len(args)==0:
		print "This command need the number of the graph that you want save and a name for the output file."
		return
	
	a=args.split()
	
	whichplot=int(a[0])
	whatset=int(a[1])
	outfile=a[2]
	plot=self._get_displayed_plot()
	#print plot
	
        print 'Select two points'
        points=self._measure_N_points(N=2, whatset=whatset)
	minbound=min(points[0].index, points[1].index)
	maxbound=max(points[0].index, points[1].index)

	xarr=[]
	yarr=[]
	try:
	  dataset=self.plots[whichplot].vectors[whatset]
	except:
          print "Invalid whichplot."
          return
        
	xarr=dataset[0][minbound:maxbound]
	yarr=dataset[1][minbound:maxbound]



	f=open(outfile,'w+')
	for i in range(len(yarr)):
		f.write(str(xarr[i])+";"+str(yarr[i])+"\n")
        f.close()

