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
    _,surface_idx,_ = Analysis.get_surface_index(ex.retract,n_smooth=300,
                                                 last_less_than=False)
    # use hooke on the data...
    fitter = flatfilts.flatfiltsCommands()
    convfilt_config = dict(convolution=[1,1,1,0,0,0],
                           blindwindow=1e3,
                           stable=0.005,
                           positive=True,
                           maxcut=0.2,
                           seedouble=300,
                           mindeviation=0)
    fitter.convfilt_config = convfilt_config
    fitter.find_contact_point = lambda *args : surface_idx
    peaks,peaks_size = fitter.has_peaks(hooke_object(time,list(force)))
    # convert force to pN for this example
    force *= 1e12
    idx_events = ex.get_retract_event_idx()
    time_events = [time[s] for s in idx_events]
    force_events = [force[s] for s in idx_events]
    # set up an array where the events are
    events,events_predicted = np.zeros(time.size),np.zeros(time.size)
    for s_true in idx_events:
        events[s] = 1
    for s_predicted in peaks:
        events_predicted[s_predicted] = 1
    fig = PlotUtilities.figure()
    plt.subplot(2,1,1)
    plt.plot(time,force,color='k',alpha=0.3)
    for t,f in zip(time_events,force_events):
        plt.plot(t,f,color='r')
    for p in peaks:
        plt.plot(time[p],force[p],'b.')
    PlotUtilities.lazyLabel("","Force(pN)","")
    plt.subplot(2,1,2)
    plt.plot(time,events,linewidth=4,label="True events",color='b')
    plt.plot(time,events_predicted,label="Predicted",color='r')
    PlotUtilities.lazyLabel("Time (s)","Force(pN)","")
    PlotUtilities.savefig(fig,"./out_hooke.png")

if __name__ == "__main__":
    run()
