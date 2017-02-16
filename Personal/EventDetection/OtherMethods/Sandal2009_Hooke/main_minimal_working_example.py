# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from Research.Personal.EventDetection.OtherMethods import method_helper
from Research.Personal.EventDetection.Util import Analysis

from GeneralUtil.python import PlotUtilities
from Research.Personal.EventDetection.OtherMethods.Sandal2009_Hooke.\
    hooke_src.trunk import flatfilts

class hooke_object:

    def __init__(self,time,force):
        self.vectors = [ [],[time,force]]

def call_hooke(split_fec,convolution=[1,1,1,0,0,0],blindwindow=None,
               seedouble=None,stable=None,positive=True,maxcut=0.5,
               mindeviation=0.5):
    """
    Args:
        split_fec: instance of split_force_extension to use
    
        Convolution:the shape of the peak to use to convolve with our vector,
        used in 'libpeakspot.conv_dx'
        
        blindwindow: the amount of x*1e9, relative to the start, to exclude. 
        essentially, another mask. It must be given in units of (1e9 * x)
        (I assume this is because 'flatfilts.has_peaks' thinks x is in nm, and 
        blindwindow is in meters). If none, defaults to a single autocorrelation
        time 
        
        stable: used in 'libpeakspot.noise_absdev', abs deviation thresh:
        'we cut the most negative (or positive) data points until the absolute 
        deviation becomes stable (it doesn't vary more than stable) or we have 
        cut more than maxcut*len(data) points.' If None, defaults to using
        the standard deviation from a correlation-time smoothed spline 
        interpolation
        
        positive: used in 'libpeakspot.noise_absdev'; by default, we assume 
        positive points are interesting (ie: force on the molecule, with + 
        being pulled away from the surface). Default is true, data has been
        zeroed
        
        maxcut: used in 'libpeakspot.noise_absdev', the maximum fraction of
        points we mask out as being below the noise floor:
        'we cut the most negative (or positive) data points until the absolute 
        deviation becomes stable (it doesn't vary more than stbale) or we have 
        cut more than maxcut*len(data) points.'
        
        seedouble: in 'libpeakspot.find_peaks' how we remove double peaks; the 
        min number of points between adjacent peaks to be considered distinct. 
        'value at which we want to "delete" double peaks. That is, if two peaks 
        have a distance < than $seedouble points , only the first is kept.'
        
        mindeviation: passed to libpeakspot.abovenoise, which 
        'Generates a vector which is 0 where the vector is less than 
        abs_devs*noise_level ; 1 if not (spike).'
        In other words, this is a noise floor mask
        
    Returns:
        see fitter.has_peaks
    """
    time,force = split_fec.retract.Time, split_fec.retract.Force                                                 
    # determine the parameters we need
    tau_time = split_fec.tau
    tau_num_points = split_fec.tau_num_points
    interpolator = Analysis.spline_interpolator(tau_time,time,force)
    force_interpolated = interpolator(time)
    mean,stdev = Analysis.spline_residual_mean_and_stdev(force,force_interpolated)
    # determine some of the parameters from what we care about
    if (blindwindow is None):
        blindwindow = tau_time 
    if (seedouble is None):
        seedouble = tau_num_points
    if (stable is None):
        stable=0.005
    convfilt_config = dict(convolution=convolution,
                           blindwindow=blindwindow,
                           stable=stable,
                           positive=positive,
                           maxcut=maxcut,
                           mindeviation=mindeviation,
                           seedouble=seedouble)
    print(convfilt_config)
    fitter = flatfilts.flatfiltsCommands()                           
    fitter.convfilt_config = convfilt_config
    _,surface_idx,_ = Analysis.get_surface_index(split_fec.retract,
                                                 n_smooth=tau_num_points,
                                                 last_less_than=not positive)
    fitter.find_contact_point = lambda *args : surface_idx
    peaks_predicted,peaks_size = fitter.has_peaks(hooke_object(time,list(force)))
    return peaks_predicted,peaks_size
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    ex = method_helper.get_example()
    retract = ex.retract
    time,force = retract.Time,retract.Force
    peaks_predicted,peaks_size = call_hooke(ex)
    print(peaks_predicted)
    # convert force to pN for this example
    force *= 1e12
    idx_events = ex.get_retract_event_idx()
    time_events = [time[s] for s in idx_events]
    force_events = [force[s] for s in idx_events]
    # set up an array where the events are
    events,events_predicted = np.zeros(time.size),np.zeros(time.size)
    for s_true in idx_events:
        events[s] = 1
    for s_predicted in peaks_predicted:
        events_predicted[s_predicted] = 1
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    plt.plot(time,force,color='k',alpha=0.3)
    for t,f in zip(time_events,force_events):
        plt.plot(t,f,color='r')
    for p in peaks_predicted:
        plt.plot(time[p],force[p],'b.')
    PlotUtilities.lazyLabel("","Force(pN)","")
    plt.subplot(2,1,2)
    plt.plot(time,events,linewidth=4,label="True events",color='b')
    plt.plot(time,events_predicted,label="Predicted",color='r',alpha=0.3)
    PlotUtilities.lazyLabel("Time (s)","Force(pN)","")
    PlotUtilities.savefig(fig,"./out_hooke.png")

if __name__ == "__main__":
    run()
