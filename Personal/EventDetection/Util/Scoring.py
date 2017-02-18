# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

    
from sklearn import metrics

class score:
    def __init__(self,x,idx_true,idx_predicted):
        """
        gets the scores associated with the binary 'is there an event happening'
        matrices given by true and predicted, which index into x
        
        Args:
            x: the data values associated with the events we care about 
            true: matrix of same length as x, 0/1 where an event is/isnt tagged
            predicted: like true, but for prediction
        """
        # save where the events are
        self.idx_true = idx_true
        self.idx_predicted = idx_predicted
        # get the 0/1 booleans arrays
        true,predicted = self.get_boolean_arrays(x)
        # have x start at 0...
        self.precision = metrics.precision_score(true,predicted)
        self.recall = metrics.precision_score(true,predicted)
        # get the x values where we think events are happenings
        true_x = x[np.where(true)]
        predicted_x = x[np.where(predicted)]
        # get the minimum distance for each
        closest_true = lambda x: true_x[np.argmin(np.abs(true_x-x))]
        self.minimum_distance_distribution = [np.abs(x-closest_true(x))
                                              for x in predicted_x]
        self.minimum_distance_median = \
            np.median(self.minimum_distance_distribution)
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
    idx_events = split_fec_with_events.get_retract_event_centers()
    separation = split_fec_with_events.retract.Separation
    return score(separation,idx_events,idx_predicted_centers)
