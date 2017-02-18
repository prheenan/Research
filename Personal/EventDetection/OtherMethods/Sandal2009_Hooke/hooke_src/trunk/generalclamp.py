#!/usr/bin/env python

'''
GENERALCLAMP.py

Plugin regarding general force clamp measurements
'''
from libhooke import WX_GOOD, ClickedPoint
import wxversion
import libhookecurve as lhc
wxversion.select(WX_GOOD)
from wx import PostEvent

class generalclampCommands:

    def plotmanip_clamp(self, plot, current, customvalue=False):
        '''
        Handles some viewing options for the "force clamp" data format, depending on the state of these configuration variables:
        (1) If self.config['fc_showphase'] != 0, the 'phase' data column (i.e. the 2nd) is shown in the 0th graph (else it isn't)
        (2) If self.config['fc_showimposed'] != 0, the 'imposed deflection' data column (i.e. the 5th) is shown in the 1st graph (else it isn't)
        (3) If self.config['fc_interesting'] == 0, the entire curve is shown in the graphs; if it has a non-zero value N, only phase N is shown.

        NOTE - my implementation of point(3) feels quite awkward - someone smarter than me plz polish that!

        '''
        
        #not a fclamp curve...
        if current.curve.experiment != 'clamp':
            return plot

        if self.config['fc_interesting'] != 0 and plot.destination==0:
            lower = int((self.config['fc_interesting'])-1)
            upper = int((self.config['fc_interesting'])+1)
            trim = current.curve.trimindexes()[lower:upper]
            newtime = []
            newzpiezo = []
            newphase = []
            for x in range(trim[0],trim[1]):
                newtime.append(self.plots[0].vectors[0][0][x])
                newzpiezo.append(self.plots[0].vectors[0][1][x])
                newphase.append(self.plots[0].vectors[1][1][x])
            self.plots[0].vectors[0][0] = newtime
            self.plots[0].vectors[0][1] = newzpiezo
            self.plots[0].vectors[1][0] = newtime
            self.plots[0].vectors[1][1] = newphase

        if self.config['fc_interesting'] != 0 and plot.destination==1:
            lower = int((self.config['fc_interesting'])-1)
            upper = int((self.config['fc_interesting'])+1)
            trim = current.curve.trimindexes()[lower:upper]
            newtime = []
            newdefl = []
            newimposed = []
            for x in range(trim[0],trim[1]):
                newtime.append(self.plots[1].vectors[0][0][x])
                newdefl.append(self.plots[1].vectors[0][1][x])
                newimposed.append(self.plots[1].vectors[1][1][x])
            self.plots[1].vectors[0][0] = newtime
            self.plots[1].vectors[0][1] = newdefl
            self.plots[1].vectors[1][0] = newtime
            self.plots[1].vectors[1][1] = newimposed            
                        
        if self.config['fc_showphase'] == 0 and plot.destination==0:
            self.plots[0].remove_set(1)
            
        if self.config['fc_showimposed'] == 0 and plot.destination==1:
            self.plots[1].remove_set(1)
                         
        return plot
      
    def do_time(self,args):
        '''
        Measures the time difference (in seconds) between two points
        Implemented only for force clamp
        ----
        Syntax: time
        '''
        if self.current.curve.experiment == 'clamp':
            time=self._delta(set=0)[0]
            print str(time*1000)+' ms'
        else:
            print 'This command makes no sense for a non-force clamp experiment.'
            
    def do_zpiezo(self,args):
        '''
        Measures the zpiezo difference (in nm) between two points
        Implemented only for force clamp
        ----
        Syntax: zpiezo
        '''
        if self.current.curve.experiment == 'clamp':
            zpiezo=self._delta(set=0)[2]
            print str(zpiezo*(10**9))+' nm'
        else:
            print 'This command makes no sense for a non-force clamp experiment.'
            
    def do_defl(self,args):
        '''
        Measures the deflection difference (in nm) between two points
        Implemented only for force clamp
        NOTE: It makes sense only on the time VS defl plot; it is still not masked for the other plot...
        -----
        Syntax: defl
        '''
        if self.current.curve.experiment == 'clamp':
            print "Warning - don't use on the zpiezo plot!"
            defl=self._delta(set=1)[2]
            print str(defl*(10**12))+' pN'
        else:
            print 'This command makes no sense for a non-force clamp experiment.'
            
    def do_step(self,args):
        '''
        Measures the length and time duration of a time-Z step
        -----
        Syntax: step
        '''
        if self.current.curve.experiment == 'clamp':
            print 'Click three points in this fashion:'
            print ' (0)-------(1)'
            print '           |'
            print '           |'
            print '           (2)----------'
            points=self._measure_N_points(N=3,whatset=0)
            dz=abs(points[2].graph_coords[1]-points[1].graph_coords[1])*(10e+8)
            dt=abs(points[1].graph_coords[0]-points[0].graph_coords[0])
            print 'dZ: ',dz,' nm'
            print 'dT: ',dt,' s'
            
        else:
            print 'This command makes no sense for a non-force clamp experiment.'

    def do_fcfilt(self,args):
        '''
        Filters out featureless force clamp curves of the current playlist.
        It's very similar to 'flatfilt' for velocity clamp curves.
        Creates a new playlist only containing non-empty curves.

        WARNING - Only works if you set an appropriate fc_interesting config variable!
        WARNING - arguments are NOT optional at the moment!

        Syntax: fcfilt maxretraction(nm) mindeviation (pN)

        Suggested values for an (i27)8 experiment with our setup are 200nm and 10-15 pN
        '''

        if self.config['fc_interesting'] == 0:
            print 'You must specify the phase of interest (using set fc_interesing X) prior to running fcfilt!'
            return
        
        maxretraction=0
        threshold=0
        args=args.split(' ')
        if len(args)==2:
            maxretraction=int(args[0])
            threshold=int(args[1])
        else:
            print 'Arguments are not optional for fcfilt. You should pass two numbers:'
            print '(1) the maximum plausible piezo retraction in NANOMETERS (e.g. the length of the protein)'
            print "(2) the threshold, in PICONEWTONS. If signal deviates from imposed more than this, it's an event"
            return
        

        print 'Processing playlist... go get yourself a cup of coffee.'
        notflat_list=[]

        c=0

        for item in self.current_list:
            c+=1
            try:
                notflat=self.has_stuff(item,maxretraction,threshold)
                print 'Curve',item.path,'is',c,'of',len(self.current_list),'--->Has Stuff =',notflat
            except:
                notflat=False
                print 'Curve',item.path,'is',c,'of',len(self.current_list),'--->could not be processed'
            if notflat:
                item.features=notflat
                item.curve=None
                notflat_list.append(item)

        if len(notflat_list)==0:
            print 'Nothing interesting here. Reconsider either your filtering criteria or your experimental data'
            return
        else:
            print 'Found',len(notflat_list),'potentially interesting curves.'
            print 'Regenerating Playlist...'
            self.pointer=0
            self.current_list=notflat_list
            self.current=self.current_list[self.pointer]
            self.do_plot(0)

    def has_stuff(self,item,maxretraction,threshold):
        '''
        Decides whether a curve has some features in the interesting phase.
        Algorithm:
            - clip the interesting phase portion of the curve.
            - discard the first 20 milliseconds (this is due to a quirk of our hardware).
            - look at the zpiezo plot and note down when (if) retratcs more than [maxretraction] nm away from the first point.
            - clip off any data after this point, with an excess of 100 points (again, an hardware quirk)
            - if the remainder is less than 100 points, ditch the curve.
            - now look at the deflection plot and check if there are points more than [threshold] pN over the 'flat zone'.
            - if you find such points, bingo!            
        '''

        item.identify(self.drivers)
   
        lower = int((self.config['fc_interesting'])-1)
        upper = int((self.config['fc_interesting'])+1)
        trim_idxs = item.curve.trimindexes()[lower:upper]
        lo=trim_idxs[0]+20                                                  #clipping the first 20 points off...
        hi=trim_idxs[1]
        trimmed_zpiezo=item.curve.default_plots()[0].vectors[0][1][lo:hi]
        trimmed_defl=item.curve.default_plots()[1].vectors[0][1][lo:hi]
        trimmed_imposed=item.curve.default_plots()[1].vectors[1][1][lo:hi]
        imposed=trimmed_imposed[21]                                         #just to match the 20-pts clipping...
        
        item.curve.close_all()
        del item.curve
        del item

        starting_z=trimmed_zpiezo[0]
        plausible=starting_z-(maxretraction*1e-9)
        det_trim=0
        while trimmed_zpiezo[det_trim]>plausible:
            det_trim+=1
            if det_trim >= len(trimmed_zpiezo):                              #breaking cycles makes me shiver...
                det_trim=len(trimmed_zpiezo)                                 #but I cannot think of anything better now.
                break
        further_trim=det_trim-100
        if further_trim<100:
            return False
        trimmed_defl=trimmed_defl[:further_trim]

        trimmed_defl.sort()
        ninetypercent=int(0.9*len(trimmed_defl))
        j=0
        sum=0
        for j in trimmed_defl[:ninetypercent]:
            sum+=j
        avg=float(sum/ninetypercent)
        sweetspot=float(avg+(threshold*1e-12))
        if trimmed_defl[-1]>sweetspot:
            flag=True
        else:
            flag=False

        del trimmed_defl,trimmed_zpiezo,trimmed_imposed            

        return flag            
        