# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

    
from sklearn import metrics
from Research.Personal.EventDetection.Util import Analysis

def get_true_and_predicted_rupture_information(example_split,
                                               idx_true,idx_predicted,
                                               num_tau=2):
    """
    Returns the rupture information (force and loading rate) at n_points
    before the indices given by this score

    Args:
        time/force: from the force extension curve
        idx_<true/predicted>: where the true/predictd events are 
        happening, zeroed to index zero of the *retract*
        num_tau: number of auto correlation times before the events to fit
    Returns:
        tuple: <list of true, list of predicted> rupture objects
    """
    retract = example_split.retract
    time,force = retract.Time,retract.Force
    n = time.size
    n_points = 2*example_split.tau_num_points
    m_slice = lambda event_idx: slice(max(event_idx-n_points,0),
                                      event_idx,1)
    rupture_func = lambda slice_ev: \
        Analysis.loading_rate_rupture_force_and_index(time,force,slice_ev)
    # need at least two points to fit a line for the predicted values
    condition = lambda slice_ev: time[slice_ev].size > 1
    true = [rupture(*rupture_func(m_slice(e))) for e in idx_true]
    pred = [rupture(*rupture_func(m_slice(e))) 
            for e in idx_predicted if condition(m_slice(e))]
    return true,pred

class rupture:
    def __init__(self,loading_rate,rupture_force,index):
        self.loading_rate = loading_rate
        self.rupture_force = rupture_force
        self.index = index

class score:
    def __init__(self,split_fec,idx_predicted):
        """
        gets the scores associated with the binary 'is there an event happening'
        matrices given by true and predicted, which index into x
        
        Args:
            split_fec: the split force extension curve with the data we 
            care about. 
            idx_predicted: where we *think* an evetnt is happening
        Returns:
            a scoring object; informaiton on precision, recall relative to
            ground truth, as well as loading rate and rupture force
        """
        idx_true = split_fec.get_retract_event_centers()
        retract = split_fec.retract
        self.tau_num_points = split_fec.tau_num_points
        x = retract.Separation
        # save where the events are
        self.idx_true = idx_true
        self.idx_predicted = idx_predicted
        # how many events are there?
        self.n_events = len(idx_true)
        self.n_events_predicted = len(idx_predicted)
        try:
            self.true_x = [x[i] for i in idx_true]
            self.pred_x = [x[i] for i in idx_predicted]
            raise IndexError("foobar")
        except IndexError as e:
            print("FEC {:s} was messed up, events at (absolute): {:s}".\
                   format(retract.Meta.Name,str(split_fec.retract.Events))
        # XXX f that
        self.ruptures_true,self.ruptures_predicted = \
            get_true_and_predicted_rupture_information(split_fec,
                                                       self.idx_true,
                                                       self.idx_predicted)
    def get_boolean_arrays(self,time):
        """
        Returns:
            arrays of the same length as the x values, 1 where there is an 
            event, 0 where there isnt. <true,predicted>
        """
        events,events_predicted = np.zeros(time.size),np.zeros(time.size)
        for s_true in self.idx_true:
            events[s_true] = 1
        for s_predicted in self.idx_predicted:
            events_predicted[s_predicted] = 1
        return events,events_predicted

        
        
        
def get_scoring_info(split_fec_with_events,idx_predicted_centers):
    return score(split_fec_with_events,idx_predicted_centers)
