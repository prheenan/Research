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

def m_colorbar(im,fig):
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    PlotUtilities.colorbar(label="Height (nm)",
                           bar_kwargs=dict(mappable=im,cax=cbar_ax))


def prep_figure(src_gel_photo,src_dna_image):
    # read in the ibw wave
    wave = CheckpointUtilities.getCheckpoint("image.pkl",
                                             PxpLoader.read_ibw_as_wave,False,
                                             src_dna_image)
    qlow,qhigh = np.percentile(wave.height_nm_rel(),[30,99])
    # range will be in nanomeyers
    range_plot = wave.range_microns() * 1000
    image = CheckpointUtilities.getCheckpoint("gel.pkl",
                                              plt.imread,False,
                                              src_gel_photo)
    plt.subplot(1,2,1)
    plt.imshow(image,extent=[0,1,0,1])
    PlotUtilities.FormatImageAxis()
    # add an arrow at the 2KB line
    ax = plt.gca()
    ax.annotate('2kbp standard', fontsize=20, xy=(.11, .57),
                xycoords='data', xytext=(50, -150),
                textcoords='offset points',
                arrowprops=dict(width = 5.,
                                headwidth = 20.,
                                headlength=20,
                                shrink = 0.05,
                                linewidth = 2,
                                alpha=0.5,
                                color = 'orange'),
                bbox=dict(boxstyle="round", fc="orange",alpha=0.3)
            )
    PlotUtilities.lazyLabel("","","Purification of 647 nm DNA")
    plt.subplot(1,2,2)
    kwargs = dict(vmin=qlow,vmax=qhigh,cmap=plt.cm.afmhot,range_plot=range_plot)
    im = ImageUtil.PlotImage(wave,aspect='equal',**kwargs)
    PlotUtilities.lazyLabel("x position (nm)","y position (nm)",
                            "AFM image of mica-bound DNA")
    ax = plt.gca()
    ax.invert_yaxis()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right') 
    return im


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
    subplots_adjust = dict(left=0.07,right=0.7,hspace=0.0,bottom=0.0,top=1.0)
    # make the figure for the presentaiton, without subplot labels and larger
    fig = PlotUtilities.figure((16,10))
    im = prep_figure(src_gel_photo,src_dna_image)
    plt.subplots_adjust(**subplots_adjust)
    m_colorbar(im,fig)
    PlotUtilities.savefig(fig,"./prep_pres.pdf",
                          subplots_adjust=subplots_adjust)
    # make the figure for the paper, with subplot labels
    fig = PlotUtilities.figure((16,8))
    im = prep_figure(src_gel_photo,src_dna_image)
    plt.tight_layout()
    m_colorbar(im,fig)
    axis_func = lambda x: x[:-1]
    PlotUtilities.label_tom(fig,loc=(0,1.1),axis_func=axis_func)
    PlotUtilities.savefig(fig,"./prep.pdf",
                          subplots_adjust=subplots_adjust)


if __name__ == "__main__":
    run()
