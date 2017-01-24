# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.EnergyLandscapes import IWT_Util,IWT_Plot
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from IgorUtil.PythonAdapter import PxpLoader
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass
from GeneralUtil.python import CheckpointUtilities,PlotUtilities

def split(RawData,idx_end_of_unfolding,fraction_for_vel=0.5):
    unfold,refold = [],[]
    for d in RawData:
        kwargs = dict(idx_end_of_unfolding=idx_end_of_unfolding,
                      fraction_for_vel=fraction_for_vel,flip_forces=True)
        unfold_tmp,refold_tmp = IWT_Util.split_into_iwt_objects(d,**kwargs)
        unfold.append(unfold_tmp)
        refold.append(refold_tmp)
    return unfold,refold

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    Base = "./"
    OutBase = Base + "out/"
    InFiles = [Base + "PatrickIsGreedy.pxp"]
    RawData = IWT_Util.\
              ReadInAllFiles(InFiles,Limit=50,
                             ValidFunc=PxpLoader.valid_fec_allow_endings)
    # get the start/ends of the re-folding and unfolding portions
    # hard coded constant for now...
    # XXX for re-folding, need to add in schedule
    # XXX ensure properly zeroed?
    idx_end_of_unfolding = int(16100/2)
    IwtData,IwtData_fold = split(RawData,idx_end_of_unfolding)
    # get the titled landscape...
    all_landscape = [-np.inf,np.inf]
    Bounds = IWT_Util.BoundsObj(bounds_folded_nm= all_landscape,
                                bounds_transition_nm= all_landscape,
                                bounds_unfolded_nm=all_landscape,
                                force_one_half_N=15e-12)
    OutBase = "./out/"
    # get the unfolding histograms
    forces_unfold = np.concatenate([r.Force for r in IwtData])
    separations_unfold = np.concatenate([r.Extension for r in IwtData])
    # get the folding histograms..
    forces_fold = np.concatenate([r.Force for r in IwtData_fold])
    separations_fold = np.concatenate([r.Extension for r in IwtData_fold])
    # zero everything...
    n_bins = 80
    fig = PlotUtilities.figure(figsize=(12,16))
    kwargs_histogram = dict(AddAverage=False,nBins=n_bins)
    plt.subplot(2,1,1)
    IWT_Util.ForceExtensionHistograms(separations_unfold*1e9,
                                      forces_unfold*1e12,**kwargs_histogram)
    PlotUtilities.xlabel("")
    PlotUtilities.title("*Unfolding* 2-D histogram")
    plt.subplot(2,1,2)
    IWT_Util.ForceExtensionHistograms(separations_fold*1e9,
                                      forces_fold*1e12,**kwargs_histogram)
    PlotUtilities.title("*Folding* 2-D histogram")
    PlotUtilities.savefig(fig,OutBase + "0_{:d}hist.pdf".format(n_bins))
    IWT_Plot.InTheWeedsPlot(OutBase=OutBase,
                            UnfoldObj=IwtData,
                            bounds=Bounds,Bins=[40,60,80,120,200],
                            max_landscape_kT=None,
                            min_landscape_kT=None)


if __name__ == "__main__":
    run()
