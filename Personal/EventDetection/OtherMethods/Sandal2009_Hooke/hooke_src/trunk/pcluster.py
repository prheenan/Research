#!/usr/bin/env python
# -*- coding: utf-8 -*-

from mdp import pca
from libhooke import WX_GOOD, ClickedPoint
import wxversion
wxversion.select(WX_GOOD)
from wx import PostEvent
import numpy as np
import scipy as sp
import copy
import os.path
import time
import libhookecurve as lhc
import pylab as pyl

import warnings
warnings.simplefilter('ignore',np.RankWarning)


class pclusterCommands:
    
    def _plug_init(self):
        self.clustplot1=None
        self.clustplot2=None

    def do_pcluster(self,args):
        '''
        pCLUSTER
        (pcluster.py)
        Automatically measures peaks and extracts informations for further clustering
        (c)Paolo Pancaldi, Massimo Sandal 2009
        '''
        if self.config['hookedir'][0]=='/':
            slash='/' #a Unix or Unix-like system
        else:
            slash='\\'
        blindw = str(self.convfilt_config['blindwindow'])
        pclus_dir = "pCluster_blind"+blindw+"_"+time.strftime("%Y%m%d_%H%M")
        self.my_work_dir = os.getcwd()+slash+pclus_dir+slash
        self.my_curr_dir = os.path.basename(os.getcwd())
        os.mkdir(self.my_work_dir)
        
        #--Custom persistent length
        pl_value=None
        for arg in args.split():
            #look for a persistent length argument.
            if 'pl=' in arg:
                pl_expression=arg.split('=')
                pl_value=float(pl_expression[1]) #actual value
            else:
                pl_value=None
                
        #configuration variables
        min_npks = self.convfilt_config['minpeaks']
        min_deviation = self.convfilt_config['mindeviation']
        
        pclust_filename = "automeasure_"+self.my_curr_dir+"_blind"+blindw+".txt" #raw_input('Automeasure filename? ')
        realclust_filename = "coordinate_"+self.my_curr_dir+"_blind"+blindw+".txt" #raw_input('Coordinates filename? ')
        peackforce_filename = "peakforce_"+self.my_curr_dir+"_blind"+blindw+".txt"  #raw_input('Peacks and Forces filename? ')
        
        f=open(self.my_work_dir+pclust_filename,'w+')
        f.write('Analysis started '+time.asctime()+'\n')
        f.write('----------------------------------------\n')
        f.write('; Contour length (nm)  ;  Persistence length (nm) ;  Max.Force (pN)  ;  Slope (N/m) ;  Sigma contour (nm) ; Sigma persistence (nm)\n')
        f.close()
        
        f=open(self.my_work_dir+realclust_filename,'w+')
        f.write('Analysis started '+time.asctime()+'\n')
        f.write('----------------------------------------\n')
        f.write('; Peak number ; Mean delta (nm)  ;  Median delta (nm) ;  Mean force (pN)  ;  Median force (pN) ; First peak length (nm) ; Last peak length (nm) ; Max force (pN) ; Min force (pN) ; Max delta (nm) ; Min delta (nm) ; Peaks Diff\n')
        f.close()
        
        f=open(self.my_work_dir+peackforce_filename,'w+')
        f.write('Analysis started '+time.asctime()+'\n')
        f.write('----------------------------------------\n')
        f.write('; Peak number  ;  1 peak Length (nm) ; 1 peak Force (pN) ;  2 peak Length (nm) ; 2 peak Force (pN) ;  3 peak Length (nm) ; 3 peak Force (pN) ;  4 peak Length (nm) ; 4 peak Force (pN) ;  5 peak Length (nm) ; 5 peak Force (pN) ;  6 peak Length (nm) ; 6 peak Force (pN) ;  7 peak Length (nm) ; 7 peak Force (pN) ;  8 peak Length (nm) ; 8 peak Force (pN)\n')
        f.close()
        

        def plot_informations(itplot,pl_value):
            '''
            OUR VARIABLES
            contact_point.absolute_coords		(2.4584142802103689e-007, -6.9647135616234017e-009)
            peak_point.absolute_coords			(3.6047748250571423e-008, -7.7142802788854212e-009)
            other_fit_point.absolute_coords	(4.1666139243838867e-008, -7.3759393477579707e-009)
            peak_location										[510, 610, 703, 810, 915, 1103]
            peak_size												[-1.2729111505202212e-009, -9.1632775347399312e-010, -8.1707438353929907e-010, -8.0335812578148904e-010, -8.7483955226387558e-010, -3.6269619757067322e-009]
            params													[2.2433999931959462e-007, 3.3230248825175678e-010]
            fit_errors											[6.5817195369767644e-010, 2.4415923138871498e-011]
            '''
            fit_points=int(self.config['auto_fit_points']) # number of points to fit before the peak maximum <50>
            
            T=self.config['temperature'] #temperature of the system in kelvins. By default it is 293 K. <301.0>
            cindex=self.find_contact_point(itplot[0]) #Automatically find contact point <158, libhooke.ClickedPoint>
            contact_point=self._clickize(itplot[0].vectors[1][0], itplot[0].vectors[1][1], cindex)
            self.basepoints=[]
            base_index_0=peak_location[-1]+self.fit_interval_nm(peak_location[-1], itplot[0], self.config['auto_right_baseline'],False)
            self.basepoints.append(self._clickize(itplot[0].vectors[1][0],itplot[0].vectors[1][1],base_index_0))
            base_index_1=self.basepoints[0].index+self.fit_interval_nm(self.basepoints[0].index, itplot[0], self.config['auto_left_baseline'],False)
            self.basepoints.append(self._clickize(itplot[0].vectors[1][0],itplot[0].vectors[1][1],base_index_1))
            self.basecurrent=self.current.path
            boundaries=[self.basepoints[0].index, self.basepoints[1].index]
            boundaries.sort()
            to_average=itplot[0].vectors[1][1][boundaries[0]:boundaries[1]] #y points to average
            avg=np.mean(to_average)
            return fit_points, contact_point, pl_value, T, cindex, avg

        def features_peaks(itplot, peak, fit_points, contact_point, pl_value, T, cindex, avg):
            '''
            calculate informations for each peak and add they in 
            c_lengths, p_lengths, sigma_c_lengths, sigma_p_lengths, forces, slopes
            '''
            c_leng=None
            p_leng=None
            sigma_c_leng=None
            sigma_p_leng=None
            force=None
            slope=None
            
            delta_force=10
            slope_span=int(self.config['auto_slope_span'])
            
            peak_point=self._clickize(itplot[0].vectors[1][0],itplot[0].vectors[1][1],peak)
            other_fit_point=self._clickize(itplot[0].vectors[1][0],itplot[0].vectors[1][1],peak-fit_points)
            
            points=[contact_point, peak_point, other_fit_point]
            
            params, yfit, xfit, fit_errors = self.wlc_fit(points, itplot[0].vectors[1][0], itplot[0].vectors[1][1], pl_value, T, return_errors=True)
            
            #Measure forces
            delta_to_measure=itplot[0].vectors[1][1][peak-delta_force:peak+delta_force]
            y=min(delta_to_measure)
            #Measure slopes
            slope=self.linefit_between(peak-slope_span,peak)[0]
            #check fitted data and, if right, add peak to the measurement
            if len(params)==1: #if we did choose 1-value fit
                p_leng=pl_value
                c_leng=params[0]*(1.0e+9)
                sigma_p_leng=0
                sigma_c_leng=fit_errors[0]*(1.0e+9)
                force = abs(y-avg)*(1.0e+12)
            else: #2-value fit
                p_leng=params[1]*(1.0e+9)
                #check if persistent length makes sense. otherwise, discard peak.
                if p_leng>self.config['auto_min_p'] and p_leng<self.config['auto_max_p']:
                    '''
                    p_lengths.append(p_leng)       
                    c_lengths.append(params[0]*(1.0e+9))
                    sigma_c_lengths.append(fit_errors[0]*(1.0e+9))
                    sigma_p_lengths.append(fit_errors[1]*(1.0e+9))
                    forces.append(abs(y-avg)*(1.0e+12))
                    slopes.append(slope)     
                    '''
                    c_leng=params[0]*(1.0e+9)
                    sigma_c_leng=fit_errors[0]*(1.0e+9)
                    sigma_p_leng=fit_errors[1]*(1.0e+9)
                    force=abs(y-avg)*(1.0e+12)
                else:
                    p_leng=None
                    slope=None
            #return c_lengths, p_lengths, sigma_c_lengths, sigma_p_lengths, forces, slopes
            return  c_leng, p_leng, sigma_c_leng, sigma_p_leng, force, slope


        # ------ PROGRAM -------
        c=0
        for item in self.current_list:
            c+=1
            item.identify(self.drivers)
            itplot=item.curve.default_plots()
            flatten=self._find_plotmanip('flatten') #extract flatten plot manipulator
            itplot[0]=flatten(itplot[0], item, customvalue=1)
            try:
                peak_location,peak_size=self.exec_has_peaks(item,min_deviation)
            except: 
                #We have troubles with exec_has_peaks (bad curve, whatever).
                #Print info and go to next cycle.
                print 'Cannot process ',item.path
                continue 

            if len(peak_location)==0:
                print 'No peaks!'
                continue

            fit_points, contact_point, pl_value, T, cindex, avg = plot_informations(itplot,pl_value)
            
            print '\n\nCurve',item.path, 'is',c,'of',len(self.current_list),': found '+str(len(peak_location))+' peaks.'

            #initialize output data vectors
            c_lengths=[]
            p_lengths=[]
            sigma_c_lengths=[]
            sigma_p_lengths=[]
            forces=[]
            slopes=[]

            #loop each peak of my curve
            for peak in peak_location:
                c_leng, p_leng, sigma_c_leng, sigma_p_leng, force, slope = features_peaks(itplot, peak, fit_points, contact_point, pl_value, T, cindex, avg)
                for var, vector in zip([c_leng, p_leng, sigma_c_leng, sigma_p_leng, force, slope],[c_lengths, p_lengths, sigma_c_lengths, sigma_p_lengths, forces, slopes]):
                    if var is not None:
                        vector.append(var)

            #FIXME: We need a dictionary here...
            allvects=[c_lengths, p_lengths, sigma_c_lengths, sigma_p_lengths, forces, slopes]
            for vect in allvects:
                if len(vect)==0:
                    for i in range(len(c_lengths)):
                        vect.append(0)
            
            print 'Measurements for all peaks detected:'
            print 'contour (nm)', c_lengths
            print 'sigma contour (nm)',sigma_c_lengths
            print 'p (nm)',p_lengths
            print 'sigma p (nm)',sigma_p_lengths
            print 'forces (pN)',forces
            print 'slopes (N/m)',slopes
            
            '''
            write automeasure text file
            '''
            print 'Saving automatic measurement...'
            f=open(self.my_work_dir+pclust_filename,'a+')
            f.write(item.path+'\n')
            for i in range(len(c_lengths)):
                f.write(' ; '+str(c_lengths[i])+' ; '+str(p_lengths[i])+' ; '+str(forces[i])+' ; '+str(slopes[i])+' ; '+str(sigma_c_lengths[i])+' ; '+str(sigma_p_lengths[i])+'\n')
            f.close()
            
            peak_number=len(c_lengths)
            
            '''
            write peackforce text file
            '''
            print 'Saving automatic measurement...'
            f=open(self.my_work_dir+peackforce_filename,'a+')
            f.write(item.path+'\n')
            peackforce_info = ''
            for i in range(len(c_lengths)):
                peackforce_info = peackforce_info + ' ; ' + str(c_lengths[i]) + ' ; ' + str(forces[i])
            f.write(' ; '+str(peak_number)+peackforce_info+'\n')
            f.close()
            
            '''
            calculate clustering coordinates
            '''
            if peak_number > 1:
                deltas=[]
                for i in range(len(c_lengths)-1):
                    deltas.append(c_lengths[i+1]-c_lengths[i])
                
                delta_mean=np.mean(deltas)
                delta_median=np.median(deltas)
                
                force_mean=np.mean(forces)
                force_median=np.median(forces)
                
                first_peak_cl=c_lengths[0]
                last_peak_cl=c_lengths[-1]
                
                max_force=max(forces[:-1])
                min_force=min(forces)
                
                max_delta=max(deltas)
                min_delta=min(deltas)
                
                delta_stdev=np.std(deltas)
                forces_stdev=np.std(forces[:-1])
                    
                peaks_diff=(last_peak_cl-first_peak_cl)/peak_number
                
                print 'Coordinates'
                print 'Peaks',peak_number
                print 'Mean delta',delta_mean
                print 'Median delta',delta_median
                print 'Mean force',force_mean
                print 'Median force',force_median
                print 'First peak',first_peak_cl
                print 'Last peak',last_peak_cl
                print 'Max force',max_force
                print 'Min force',min_force
                print 'Max delta',max_delta
                print 'Min delta',min_delta
                print 'Delta stdev',delta_stdev
                print 'Forces stdev',forces_stdev
                print 'Peaks difference',peaks_diff
                
                '''
                write clustering coordinates
                '''
                f=open(self.my_work_dir+realclust_filename,'a+')
                f.write(item.path+'\n')
                f.write(' ; '+str(peak_number)+     # non considerato
                        ' ; '+str(delta_mean)+      # 0
                        ' ; '+str(delta_median)+    # 1 -
                        ' ; '+str(force_mean)+      # 2
                        ' ; '+str(force_median)+    # 3 -
                        ' ; '+str(first_peak_cl)+   # 4 -
                        ' ; '+str(last_peak_cl)+    # 5 -
                        ' ; '+str(max_force)+       # 6
                        ' ; '+str(min_force)+       # 7
                        ' ; '+str(max_delta)+       # 8
                        ' ; '+str(min_delta)+       # 9
                        ' ; '+str(delta_stdev)+     # 10
                        ' ; '+str(forces_stdev)+    # 11
                        ' ; '+str(peaks_diff)+      # 12
                        '\n')
                f.close()
                
        # start PCA
        self.do_pca(pclus_dir+"/"+realclust_filename)
        
                
    def do_pca(self,args):
        '''
        PCA -> "pca gaeta_coor_blind50.txt 1,3,6"
        Automatically measures pca from coordinates filename and shows two interactives plots
        With the second argument (arbitrary) you can select the columns and the multiplier factor 
        to use for the pca (for es "1,3*50,6,8x10,9"). Dont use spaces. "*" or "x" are the same thing.
        Without second argument it reads pca_config.txt file
        (c)Paolo Pancaldi, Massimo Sandal 2009
        '''
        
        # reads the columns of pca
        if self.config['hookedir'][0]=='/':
            slash='/' #a Unix or Unix-like system
        else:
            slash='\\'
        self.my_hooke_dir = self.config['hookedir']+slash
        #self.my_work_dir = os.getcwd()+slash+"pCluster_"+time.strftime("%Y%m%d_%H%M")+slash
        #self.my_curr_dir = os.path.basename(os.getcwd())
        conf=open(self.my_hooke_dir+"pca_config.txt")
        config = conf.readlines()
        conf.close()
        
        self.plot_myCoord = []          # tiene le coordinate prese direttamente dal file creato con pCluster
        self.plot_origCoord = []        # tiene le coordinate solo delle colonne scelte e moltiplicate per i valori scelti
        self.plot_pcaCoord = []         # tiene le due colonne della pca
        self.plot_pcaCoordTr = []       # tiene le due colonne della pca trasposta
        self.plot_FiltOrigCoord = []    # tiene le coordinate solo dei punti filtrati per densita
        self.plot_FiltPaths = []        # tiene i paths dei plot solo dei punti filtrati per densita
        self.plot_paths = []            # tiene i paths di tutti i plots
        self.plot_NewPcaCoord = []      # tiene le due colonne della pca filtrate per densita
        self.plot_NewPcaCoordTr=[]      # tiene le due colonne della pca trasposta filtrate per densita
        plot_path_temp = ""
        
        # prende in inpunt un arg (nome del file) 
        # e il secondo le colonne su cui lavorare (e' arbitrario, riceve x es "1,2,3")
        arg = args.split(" ")
        if arg[0]==args:
            file_name=args
        else:
            file_name=arg[0]
            config[0] = arg[1]
        
        # creo l'array "plot_myCoord" con tutte le coordinate dei plots
        # e l'array plot_paths con tutti i percorsi dei plots
        nPlotTot = -3 #tolgo le prime 3 righe iniziali del file
        f=open(file_name)
        rows = f.readlines()
        for row in rows:
            if row[0]!=" " and row[0]!="":
                nPlotTot = nPlotTot+1
                plot_path_temp = row
            if row[0]==" " and row.find('nan')==-1 and row.find("-1.#IND")==-1:
                row = row[row.index(";",2)+2:].split(" ; ")	# non considero la prima colonna col #picchi
                row = [float(i) for i in row]
                
                #0:Mean delta, 1:Median delta, 2:Mean force, 3:Median force, 4:First peak length, 5:Last peak length
                #6:Max delta 7:Min delta 8:Max force 9:Min force 10:Std delta 11:Std force
                if (row[0]<500 and row[1]<500 and row[2]<500 and row[3]<500 and row[4]<500 and row[5]<500 and row[6]<500 and row[7]<500 and row[8]<500 and row[9]<500 and row[10]<500 and row[11]<500):
                    if (row[0]>0 and row[1]>0 and row[2]>0 and row[3]>0 and row[4]>0 and row[5]>0 and row[6]>0 and row[7]>0 and row[8]>0 and row[9]>0 and row[10]>0 and row[11]>0):
                        #row = row[0], row[2], row[3]*3, row[6], row[7]*56, row[8]
                        self.plot_myCoord.append(row)
                        self.plot_paths.append(plot_path_temp)
        f.close()
        
        # creo l'array con alcune colonne e pure moltiplicate 
        for row in self.plot_myCoord:
            res=[]
            for cols in config[0].split(","):
                if cols.find("*")!=-1:
                    col = int(cols.split("*")[0])
                    molt = int(cols.split("*")[1])
                elif cols.find("x")!=-1:
                    col = int(cols.split("x")[0])
                    molt = int(cols.split("x")[1])
                else:
                    col = int(cols)
                    molt = 1
                res.append(row[col]*molt)
            self.plot_origCoord.append(res)
        
        # array convert, calculate PCA, transpose
        self.plot_origCoord = np.array(self.plot_origCoord,dtype='float')
        #print self.plot_origCoord.shape
        self.plot_pcaCoord = pca(self.plot_origCoord, output_dim=2)	#other way -> y = mdp.nodes.PCANode(output_dim=2)(array)
        self.plot_pcaCoordTr = np.transpose(self.plot_pcaCoord)
        pca_X=np.array(self.plot_pcaCoordTr[0],dtype='float')
        pca_Y=np.array(self.plot_pcaCoordTr[1],dtype='float')
        
        '''
        # Start section of testing with good plots                                  # 4 TESTING!
        Xsyn_1=[]
        Ysyn_1=[]        
        Xgb1_1=[]
        Ygb1_1=[]
        Xbad_1=[]
        Ybad_1=[]
        goodnamefile=open(file_name.replace("coordinate", "good"),'r')
        goodnames=goodnamefile.readlines()
        nPlotGood = len(goodnames)-2 #tolgo prima e ultima riga
        goodnames=[i.split()[0] for i in goodnames[1:]]
        
        for index in range(len(self.plot_paths)):
            if self.plot_paths[index][:-1] in goodnames:
                Xsyn_1.append(pca_X[index])
                Ysyn_1.append(pca_Y[index])
            else:
                Xbad_1.append(pca_X[index])
                Ybad_1.append(pca_Y[index])
        # Stop section of testing with good plots                                   # 4 TESTING!
        '''
        
        # print first plot
        clustplot1=lhc.PlotObject()
        clustplot1.add_set(pca_X,pca_Y)
        #clustplot1.add_set(Xbad_1,Ybad_1) # 4 TESTING!
        #clustplot1.add_set(Xsyn_1,Ysyn_1) # 4 TESTING!
        clustplot1.normalize_vectors()
        clustplot1.styles=['scatter', 'scatter','scatter']
        clustplot1.colors=[None,'red','green']
        clustplot1.destination=0
        self._send_plot([clustplot1])
        self.clustplot1=clustplot1
        
        # density and filer estimation
        kernel = sp.stats.kde.gaussian_kde(sp.c_[pca_X,pca_Y].T)
        tallest = 0
        for i in range(len(pca_X)):
            kern_value = kernel.evaluate([pca_X[i],pca_Y[i]])
            if tallest < kern_value:
                    tallest = float(kern_value)
        if float(config[1]) == 0:
            my_filter = float(tallest / 3.242311147)
        else:
            my_filter = float(config[1])
        '''
        # section useful only for graphic printing
        xmin = pca_X.min()
        xmax = pca_X.max()
        ymin = pca_Y.min()
        ymax = pca_Y.max()
        mX, mY = sp.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        Z = sp.rot90(sp.fliplr(sp.reshape(kernel(sp.c_[mX.ravel(), mY.ravel()].T).T, mX.T.shape)))
        axis_X = np.linspace(xmin,xmax,num=100)
        axis_Y = np.linspace(ymin,ymax,num=100)
        '''
        
        # density filtering:
        # tramite "kernel.evaluate" trovo lo score (altezza) di ogni coordinata e decido se mantenerla o no
        filtered_pca_X = []
        filtered_pca_Y = []
        filtered_PcaCoordTr = []
        filtered_PcaCoord = []
        for i in range(len(pca_X)):
            kern_value = kernel.evaluate([pca_X[i],pca_Y[i]])
            if kern_value > my_filter:
                filtered_pca_X.append(pca_X[i])
                filtered_pca_Y.append(pca_Y[i])
        filtered_PcaCoordTr.append(filtered_pca_X)
        filtered_PcaCoordTr.append(filtered_pca_Y)
        filtered_PcaCoord = np.transpose(filtered_PcaCoordTr)
        
        # creo i due array "plot_FiltOrigCoord" e "plot_FiltPaths" contenenti solo i dati filtrati con alta densita
        for index in range(len(self.plot_pcaCoord)):
            if self.plot_pcaCoord[index] in filtered_PcaCoord:
                self.plot_FiltOrigCoord.append(self.plot_myCoord[index])
                self.plot_FiltPaths.append(self.plot_paths[index])
        
        '''
        # START PCA#2: USELESS!!!
        
        # creo l array con alcune colonne e pure moltiplicate
        temp_coord = []
        for row in self.plot_FiltOrigCoord:
            res=[]
            for cols in config[2].split(","):
                if cols.find("*")!=-1:
                    col = int(cols.split("*")[0])
                    molt = int(cols.split("*")[1])
                elif cols.find("x")!=-1:
                    col = int(cols.split("x")[0])
                    molt = int(cols.split("x")[1])
                else:
                    col = int(cols)
                    molt = 1
                res.append(row[col]*molt)
            temp_coord.append(res)
        self.plot_FiltOrigCoord = temp_coord
                
        # ricalcolo la PCA: array convert, calculate PCA, transpose
        self.plot_FiltOrigCoord = np.array(self.plot_FiltOrigCoord,dtype='float')
        #print self.plot_FiltOrigCoord.shape
        self.plot_NewPcaCoord = pca(self.plot_FiltOrigCoord, output_dim=2)	#other way -> y = mdp.nodes.PCANode(output_dim=2)(array)
        self.plot_NewPcaCoordTr = np.transpose(self.plot_NewPcaCoord)
        pca_X2=np.array(self.plot_NewPcaCoordTr[0],dtype='float')
        pca_Y2=np.array(self.plot_NewPcaCoordTr[1],dtype='float')
        
        # Start section of testing with good plots                              # 4 TESTING!
        Xsyn_2=[]
        Ysyn_2=[]
        Xbad_2=[]
        Ybad_2=[]
        for index in range(len(self.plot_FiltPaths)):
            if self.plot_FiltPaths[index][:-1] in goodnames:
                Xsyn_2.append(pca_X2[index])
                Ysyn_2.append(pca_Y2[index])
            else:
                Xbad_2.append(pca_X2[index])
                Ybad_2.append(pca_Y2[index])
        
        # print second plot
        clustplot2=lhc.PlotObject()
        #clustplot2.add_set(pca_X2,pca_Y2)
        clustplot2.add_set(Xbad_2,Ybad_2)                                       # 4 TESTING!
        clustplot2.add_set(Xsyn_2,Ysyn_2)                                       # 4 TESTING!
        clustplot2.normalize_vectors()
        clustplot2.styles=['scatter', 'scatter','scatter']
        clustplot2.colors=[None,'red','green']
        clustplot2.destination=1
        self._send_plot([clustplot2])
        self.clustplot2=clustplot2
        '''
        
        # PRINT density plot
        clustplot2=lhc.PlotObject()
        clustplot2.add_set(filtered_pca_X,filtered_pca_Y)
        clustplot2.normalize_vectors()
        clustplot2.styles=['scatter', 'scatter','scatter']
        clustplot2.colors=[None,'red','green']
        clustplot2.destination=1
        self._send_plot([clustplot2])
        self.clustplot2=clustplot2
        
        # printing results
        config_pca1 = config[0].replace("*", "x").rstrip("\n")
        config_pca2 = config[2].replace("*", "x").rstrip("\n")
        print ""
        print "- START: "+file_name
        print "Curve totali: ", nPlotTot
        #print "Curve totali good: ", nPlotGood                                  # 4 TESTING!
        print "- FILTRO 1: 0-500 e NaN"
        print "Curve totali rimaste: ", len(self.plot_origCoord)
        #print 'Curve good rimaste: ', len(Xsyn_1)                               # 4 TESTING!
        print "- FILTRO 2: PCA:"+config_pca1+" e DENSITA:"+str(my_filter)
        print "Curve totali rimaste: ", len(self.plot_FiltOrigCoord)
        #print 'Curve good rimaste: ', len(Xsyn_2)                               # 4 TESTING!
        print "Piu alta: ", tallest
        #print "- FILTRO 3: 2'PCA:"+config_pca2
        print ""
        
        # -- exporting coordinates and plot of PCA in debug mode! --
        if config[3].find("true")!=-1:
            #1' PCA: save plot and build coordinate s file
            self.do_export(file_name.replace("coordinate_", "debug_pca1graph_").replace('.txt','_'+config_pca1) + " 0")
            f = open(file_name.replace("coordinate_", "debug_pca1coor_").replace('.txt','_'+config_pca1+'.txt'),'w')
            for i in range(len(pca_X)):
                f.write (str(i) + "\t" + str(pca_X[i]) + "\t" + str(pca_Y[i]) + "\n")
            f.close()
            #2' PCA: save plot and build coordinate s file
            #self.do_export(file_name.replace("coordinate_", "debug_pca2graph_").replace('.txt','_'+config_pca2) + " 1")
            #f = open(file_name.replace("coordinate_", "debug_pca2coor_").replace('.txt','_'+config_pca2+'.txt'),'w')
            #for i in range(len(pca_X2)):
            #    f.write (str(i) + "\t" + str(pca_X2[i]) + "\t" + str(pca_Y2[i]) + "\n")
            #f.close()
            #DENSITY: save plot
            self.do_export(file_name.replace("coordinate_", "debug_densitygraph_").replace('.txt','_'+config_pca1+'_'+str(my_filter).replace(".",",")) + " 1")
            f = open(file_name.replace("coordinate_", "debug_densitycoor_").replace('.txt','_'+config_pca1+'_'+str(my_filter).replace(".",",")+'.txt'),'w')
            for i in range(len(filtered_pca_X)):
                f.write (str(i) + "\t" + str(filtered_pca_X[i]) + "\t" + str(filtered_pca_Y[i]) + "\n")
            f.close()
            #ALL GOOD COORDINATES (without NaN and 0<x<500)
            f = open(file_name.replace("coordinate_", "debug_allgoodcoor_"),'w')
            for i in range(len(self.plot_myCoord)):
                for cel in self.plot_myCoord[i]:
                    f.write (" ; " + str(cel))
                f.write ("\n")
            f.close()
        
        # pCLUSTER SAVING!!!
        import shutil
        pcl_name = file_name.replace("coordinate_", "goodplots_").replace('.txt','_'+config_pca1+'_'+str(my_filter).replace(".",","))
        if os.path.exists(pcl_name+slash): shutil.rmtree(pcl_name)
        os.mkdir(pcl_name+slash)
        f = open(pcl_name+'.txt','w')
        for i in range(len(self.plot_FiltPaths)):
            myfile = str(self.plot_FiltPaths[i]).rstrip("\n")
            f.write (myfile+"\n")
            shutil.copy2(myfile, pcl_name)
        f.close()
        
        
    def do_multipca(self,args):
        '''
        MULTIPCA -> "multipca gaeta_coor_blind50.txt 3"
        Automatically multiply the column suggest in second argument for value between 1-100 (step of 2), 
        measures pca from coordinates filename and save the png plots.
        (c)Paolo Pancaldi, Massimo Sandal 2009
        '''
        # reads the columns of pca
        conf=open("pca_config.txt")
        config = conf.readlines() # config[0] = "1,2,3"
        conf.close()
        # cycling pca
        arg = args.split(" ")
        file_name=arg[0]
        column=str(arg[1])
        for i in range(1, 51, 1):
            self.do_pca(file_name + " " + config[0].replace(column,column+"*"+str(i),1))

    def do_doublepca(self,args):
        '''
        DOUBLEPCA -> "doublepca gaeta_coor_blind50.txt"
        Automatically it launches the pca command for all combinations with two column
        (c)Paolo Pancaldi, Massimo Sandal 2009
        '''
        # cycling pca
        arg = args.split(" ")
        file_name=arg[0]
        for i in range(1, 13):
            for j in range(1, 13):
                if i!=j:
                    self.do_pca(file_name + " " + str(i) + "," + str(j))
                    
    def do_triplepca(self,args):
        '''
        TRIPLEPCA -> "triplepca gaeta_coor_blind50.txt"
        Automatically it launches the pca command for all combinations with three column
        (c)Paolo Pancaldi, Massimo Sandal 2009
        '''
        # cycling pca
        arg = args.split(" ")
        file_name=arg[0]
        for i in range(1, 13):
            for j in range(1, 13):
                for k in range(1, 13):
                    if i!=j and i!=k and j!=k:
                        self.do_pca(file_name + " " + str(i) + "," + str(j) + "," + str(k))

    def do_pclick(self,args):
        '''
        It returns id, coordinates and file name of a clicked dot on a PCA graphic
        '''
        
        self._send_plot([self.clustplot1]) #quick workaround for BAD problems in the GUI
        print 'Click point'
        point = self._measure_N_points(N=1, whatset=0)
        indice = point[0].index
        plot_file = self.plot_paths[indice]
        dot_coord = self.plot_pcaCoord[indice]
        print "file: " + str(plot_file).rstrip()
        print "id: " + str(indice)
        print "coord: " + str(dot_coord)
        self.do_genlist(str(plot_file))
        #self.do_jump(str(plot_file))
        
    # indea iniziata e messa da parte...
    def do_peakforce(self, args):
        '''
        peackforce -> "peackforce peackforce_file.txt"
        Automatically measures peack and force plots
        (c)Paolo Pancaldi, Massimo Sandal 2009
        '''
        
        # prende in inpunt un arg (nome del file)
        file_name=args
        f=open(file_name)
        
        # scrivo un file temp
        g = open('_prove.txt','w')
        
        plot_file = '';
        rows = f.readlines()
        for row in rows:
            if row[0]=="/":
                plot_file = row
            if row[0]==" " and row.find('nan')==-1 and row.find("-1.#IND")==-1:
                # FILTRO SUI 7 PICCHI
                num_pic = int(row.split(" ; ")[1])
                if num_pic==7:
                    width_force = row.split(" ; ")
                    w1 = float(width_force[2]); f1 = float(width_force[3]);
                    w2 = float(width_force[4]); f2 = float(width_force[5]);
                    w3 = float(width_force[6]); f3 = float(width_force[7]);
                    w4 = float(width_force[8]); f4 = float(width_force[9]);
                    w5 = float(width_force[10]); f5 = float(width_force[11]);
                    w6 = float(width_force[12]); f6 = float(width_force[13]);
                    w7 = float(width_force[14]); f7 = float(width_force[15]);
                    if w1>0 and w1<1000 and w2>0 and w2<1000 and w3>0 and w3<1000 and w4>0 and w4<1000 and w5>0 and w5<1000 and w6>0 and w6<1000 and w7>0 and w7<1000:
                        score_76 = abs(32 - (w7 - w6))
                        score_65 = abs(32 - (w6 - w5))
                        score_54 = abs(32 - (w5 - w4))
                        score_43 = abs(32 - (w4 - w3))
                        score_32 = abs(32 - (w3 - w2))
                        score_21 = abs(32 - (w2 - w1))
                        writeme = str(score_76) + " --- " + str(row)
                        g.write(writeme)
        g.close()
        f.close()
        
        