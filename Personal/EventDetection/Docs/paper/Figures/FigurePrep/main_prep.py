# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys



sys.path.append("../../../../../../../")
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from IgorUtil.PythonAdapter import PxpLoader
from Research.Perkins.AnalysisUtil.Images import ImageUtil

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    base = FEC_Util.default_data_root()
    figure_base = base + "4Patrick/CuratedData/Masters_CSCI/Prep/finals/"
    src_gel_photo = figure_base + "P1010008_cropped.png"
    src_dna_image = figure_base + "2017-2-13-snl-3mM-nickel-2nM-DNA-10-mnute"+\
                    "-deposition-rinsed-with-water-then-buffer-transfer-to-"+\
                    "3mMNickel-Image0003.ibw"
    # read in the ibw wave
    wave = PxpLoader.read_ibw_as_wave(src_dna_image)
    qlow,qhigh = np.percentile(wave.height_nm_rel(),[30,99])
    # range will be in nanomeyers
    range_plot = wave.range_microns() * 1000
    image = plt.imread(src_gel_photo)
    fig = PlotUtilities.figure((16,8))
    plt.subplot(1,2,1)
    kwargs = dict(vmin=qlow,vmax=qhigh,cmap=plt.cm.afmhot,range_plot=range_plot)
    ax = ImageUtil.PlotImage(wave,**kwargs)
    PlotUtilities.lazyLabel("nanometers","nanometers",
                            "AFM image of mica-bound DNA")
    PlotUtilities.colorbar("Height (nm)")
    plt.subplot(1,2,2)
    plt.imshow(image,extent=[0,1,0,1])
    PlotUtilities.FormatImageAxis(aspect="auto")
    PlotUtilities.lazyLabel("","","Electrophoretic purification of 647nm DNA")
    PlotUtilities.savefig(fig,"./out.png")


if __name__ == "__main__":
    run()
