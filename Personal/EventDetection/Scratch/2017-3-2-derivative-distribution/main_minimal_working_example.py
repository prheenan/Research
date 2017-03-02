# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from GeneralUtil.python import GenUtilities,CheckpointUtilities
from Research.Personal.EventDetection.Util import Analysis,Plotting,Scoring
from Research.Personal.EventDetection._2SplineEventDetector import Detector

from GeneralUtil.python import PlotUtilities
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base_debug = "../../_1ReadDataToCache/debug_no_event/"
    files = GenUtilities.getAllFiles(base_debug,ext=".pkl")
    data = [CheckpointUtilities.getCheckpoint(f,None,False) for f in files]
    for fec in data:
        # get the zeroed force extension curve
        split_fec = Analysis.zero_and_split_force_extension_curve(fec)
        retract = split_fec.retract
        time = retract.Time
        interpolator_force_onto_time = split_fec.retract_spline_interpolator()
        derivative_force = interpolator_force_onto_time.derivative()(time)
        interpolated_force = interpolator_force_onto_time(time)
        min_interp,max_interp = \
            [min(interpolated_force),max(interpolated_force)]
        min_d = min(derivative_force)
        max_d = max(derivative_force)
        deriv_plot = (derivative_force-min_d)/(max_d-min_d) 
        deriv_plot = deriv_plot * (max_interp-min_interp) + min_interp
        # get the median and std of deriv
        med_deriv = np.median(deriv_plot)
        q25,q75 = np.percentile(deriv_plot,[25,75])
        iqr_region_idx = np.where( (deriv_plot <= q75) & 
                                   (deriv_plot >= q25))[0]
        std_iqr = np.std(deriv_plot[iqr_region_idx])
        probability = np.zeros(deriv_plot.size)
        probability[np.where(deriv_plot >= med_deriv - std_iqr)]  = 1
        possible_idx = np.where(deriv_plot <= med_deriv - std_iqr)
        possible_deriv = deriv_plot[possible_idx]
        k = (possible_deriv-med_deriv)/std_iqr
        probability[possible_idx]  = 1/k**2
        probability = np.minimum(probability,1)
        probability_plot = probability
        plt.subplot(2,1,1)
        plt.plot(retract.Force,color='k',alpha=0.3)
        plt.plot(interpolated_force,color='g',linewidth=4)
        plt.axhline(med_deriv,color='k')
        plt.axhline(med_deriv+std_iqr,color='b',linestyle='--')
        plt.axhline(med_deriv-std_iqr,color='b',linestyle='--')
        plt.subplot(2,1,2)
        plt.plot( probability_plot,color='r',linestyle='--')
        plt.semilogy()
        plt.show()
    


if __name__ == "__main__":
    run()
