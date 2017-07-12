# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../")
from IgorUtil.PythonAdapter import PxpLoader
from GeneralUtil.python import PlotUtilities,CheckpointUtilities,GenUtilities
from Research.Perkins.AnalysisUtil.Images import ImageUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util

import matplotlib.gridspec as gridspec
from scipy.stats import norm

from matplotlib_scalebar.scalebar import ScaleBar


    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = FEC_Util.default_data_root()
    input_directory = base + \
        "4Patrick/CuratedData/Outreach/2017-JILA-photo-contest/"
    files = GenUtilities.getAllFiles(input_directory,ext=".ibw")
    # read all the ibw files in and cache them
    images = []
    func = PxpLoader.read_ibw_as_image
    for i,f in enumerate(files):
        cache_file = "./{:d}.pkl".format(i)
        tmp = CheckpointUtilities.getCheckpoint(cache_file,func,False,f)
        images.append(tmp)
    crop = [ None,None,None,None,None,None,None,None ]
    for i,image in enumerate(images):        
        fig = PlotUtilities.figure()
        vmin,vmax = np.percentile(image.height_nm_rel(),[75,99])
        range_plot_nanometers = 1000 * image.range_microns()   
        vmin_dict = dict(vmin=vmin,vmax=vmax)
        ImageUtil.PlotImage(image,cmap=plt.cm.gray,
                            range_plot=range_plot_nanometers,**vmin_dict)
        PlotUtilities.FormatImageAxis()
        PlotUtilities.colorbar("Height (nm)")
        pixel_size_meters = image.pixel_size_meters
        pixel_size_nanometers = pixel_size_meters * 1e9
        print(pixel_size_nanometers,range_plot_nanometers)
        scalebar = ScaleBar(pixel_size_nanometers,'nm',box_alpha=0.7)
        plt.gca().add_artist(scalebar)                                    
        PlotUtilities.savefig(fig,"{:d}.pdf".format(i))
    
if __name__ == "__main__":
    run()
