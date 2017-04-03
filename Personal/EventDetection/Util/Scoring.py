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
        self.min_x = min(x)
        self.max_x = max(x)
        # save where the events are
        self.idx_true = idx_true
        self.idx_predicted = idx_predicted
        # how many events are there?
        self.n_events = len(idx_true)
        self.n_events_predicted = len(idx_predicted)
        meta = retract.Meta
        self.name = meta.Name
        self.source_file = meta.SourceFile
        try:
            self.true_x = [x[i] for i in idx_true]
            self.pred_x = [x[i] for i in idx_predicted]
            self.ruptures_true,self.ruptures_predicted = \
                get_true_and_predicted_rupture_information(split_fec,
                                                           self.idx_true,
                                                           self.idx_predicted)
        except IndexError as e:
            # we were fed a bad index in true (annotator error) or predicted
            # (a method went wacky). report it, and zero out this.
            self.true_x,self.pred_x = [],[]
            self.ruptures_true,self.ruptures_predicted = [],[]
            true_ev = str(split_fec.retract.Events)
            pred_ev = str(idx_predicted)
            n = retract.Force.size
            fec_name = "{:s}{:s}".format(meta.SourceFile,meta.Name)
            # note: the true events are absolute-offset to facilitate
            # easier finding and fixing them; predicted are relative.
            print("{:s} (N_retract={:d}), bad events (true/pred): {:s}/{:s}".\
                   format(fec_name,int(n),true_ev,pred_ev))
            print(e)
    def euclidean_rupture_spectrum_distance(self):
        spectrum_tuple = lambda x: (x.loading_rate,x.rupture_force)
        all_tuples = lambda list_v: np.array([spectrum_tuple(x) 
                                              for x in list_v])
        X = all_tuples(self.ruptures_true)
        Y = all_tuples(self.ruptures_predicted) 
        # get the distances from x to y and from y to x
        if (len(Y) == 0):
            dist_1 = []
            dist_2 = [self.max_displacement() for x in X]
        else:
            _,dist_1 = metrics.pairwise_distances_argmin_min(X=X,Y=Y)
            _,dist_2 = metrics.pairwise_distances_argmin_min(X=Y,Y=X)
        all_distances = list(dist_1) + list(dist_2)
        return all_distances
        
    def n_true_and_predicted_events(self):
        """
        returns te

        if no predicted events, this value is none

        Args:
             kwargs: passed to minimum_distance_distribution
        """
        return self.n_true(),self.n_pred()
    def n_true(self):
        return len(self.true_x)
    def n_pred(self):
        return len(self.pred_x)
    def max_displacement(self):
        to_ret = abs(self.max_x-self.min_x)
        return to_ret
    def minimum_distance_distribution(self,to_true=True,floor_is_max=False):
        """
        returns the median of the smallest distance from <predicted/true>
        to <true/predicted> if to_true is <true,false>

        if no predicted events, this value is none

        Args:
             to_true: if the distance should be from predicted *to* the true
             events (otherwise vice versa)
            
             floor_is_max: if the baseline (e.g. predicted) is empty, 
             should we return max_x for each of search.

             kwargs: passed to minimum_distance_distribution
        """
        if (to_true):
            baseline = self.true_x
            search = self.pred_x
        else:
            baseline = self.pred_x
            search = self.true_x
        if (len(baseline) == 0):
            if (floor_is_max):
                max_x = self.max_displacement()
                return [max_x for x in search]
            else:
                return []
        # POST: something in baseline
        closest_true = lambda x: baseline[np.argmin(np.abs(baseline-x))]
        min_distance_distribution = [np.abs(x-closest_true(x)) for x in search]
        return min_distance_distribution
    def minimum_distance_median(self,**kwargs):
        """
        returns the median of the smallest distance to an event

        if no predicted events, this value is none

        Args:
             kwargs: passed to minimum_distance_distribution
        """
        min_distance_distribution = self.minimum_distance_distribution(**kwargs)
        if (len(min_distance_distribution) > 0):
            return np.median(min_distance_distribution)
        else:
            return None

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
