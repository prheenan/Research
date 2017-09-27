# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,matplotlib,os,copy,scipy
from scipy.interpolate import griddata
sys.path.append("../../../../../../../../")

from GeneralUtil.python import GenUtilities,PlotUtilities,CheckpointUtilities
from Research.Perkins.AnalysisUtil.Images import PolymerTracing

def get_x_and_y_arrays(text_file,size_images_pixels):
    """
    Returns: the x and y columns (0 and 1) of text_file 
    """
    data = np.loadtxt(text_file)
    x = data[:,0]
    y = data[:,1]
    assert ((x > 1) & (y > 1)).all()
    return x,y
    
def get_contour_length(x,y):
    """
    Returns: the contour length of the line L_i=(x,y)_i
    """
    return np.sum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
    
def get_x_y_and_contour_lengths(text_files,size_images_pixels):
    """
    Gets all the worm_objects associated with text_files
    
    Args:
        text_files: list of files, each has columns like <x>,<y>
    Returns:
        list of worm_object 
    """
    to_ret = []
    for t in text_files:
        x,y = get_x_and_y_arrays(t,size_images_pixels)
        to_ret.append(PolymerTracing.worm_object(x,y,t))
    return to_ret     

def prh_hist(data,**hist_kw):
    counts,_,_ = plt.hist(data,**hist_kw)
    y = max(counts * 1.05)
    ax = plt.gca()
    try:
        plt.boxplot(data,positions=[y],vert=False,manage_xticks=False,meanline=True,
                    showmeans=True,flierprops=dict(color='k',markersize=1))
    except IndexError:
        # not enough data 
        pass
   
def print_info(dist,name):
    print("For {:s}:".format(name))
    print("\t Median: {:.1f} nm".format(np.median(dist)*1e9))
    print("\t Mean  : {:.1f} nm".format(np.mean(dist)*1e9))
    print("\t Std   : {:.1f} nm".format(np.std(dist)*1e9))


def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "in/"
    out_dir = "out/"
    ext = ".pkl"
    ext_text = ".txt"
    GenUtilities.ensureDirExists(out_dir)
    image_files = GenUtilities.getAllFiles(in_dir,ext=ext)
    text_files = GenUtilities.getAllFiles(in_dir,ext=ext_text)
    size_images_meters = 2e-6
    size_images_pixels = 512
    conversion_meters_per_px = size_images_meters / size_images_pixels  
    objs_all = []                        
    for file_name in image_files:
        file_no_ext = file_name.replace(ext,"")
        file_no_number = file_no_ext.split("_")[0]
        these_text_files= [ t for t in text_files 
                            if str(file_no_number) in str(t)]
        objs_tmp = get_x_y_and_contour_lengths(these_text_files,
                                               size_images_pixels)
        image_obj = CheckpointUtilities.lazy_load(file_name)
        im = image_obj.height_nm_rel()
        assert (im.shape[0] == im.shape[1]) 
        assert (im.shape[0] == size_images_pixels)
        # POST: dimensions are OK 
        img = PolymerTracing.tagged_image(image_obj,objs_tmp,file_no_number)
        objs_all.append(img)
    # POST: all the contour lengths are set in 'real' units ]  
    L0_protein = np.concatenate([o.L0_protein_dna() for o in objs_all])
    L0_dna = np.concatenate([o.L0_dna_only() for o in objs_all])
    print_info(L0_dna,"only DNA")
    print_info(L0_protein,"DNA+PRC2")
    n_protein = (L0_protein).size
    n_dna = L0_dna.size
    fig = PlotUtilities.figure()            
    n_str = lambda n: "\n(N={:d})".format(n)
    sanit_L0 = lambda x: x*1e6    
    L0_dna_plot = sanit_L0(L0_dna)
    L0_protein_plot = sanit_L0(L0_protein)
    xmin = np.min(np.concatenate([L0_dna_plot,L0_protein_plot]))
    xmax = np.max(np.concatenate([L0_dna_plot,L0_protein_plot]))
    n_bins = 12
    bins = np.linspace(xmin,xmax,endpoint=True,num=n_bins)
    xlim = [0,xmax]
    kw_dna = dict(color='g',alpha=0.3)
    kw_protein = dict(color='b',hatch='//',alpha=0.7)
    ax= plt.subplot(2,1,1)
    prh_hist(L0_dna_plot,normed=True,bins=bins,
             label="DNA Only" + n_str(n_dna),**kw_dna)
    PlotUtilities.lazyLabel("","P (1/microns)","")
    PlotUtilities.no_x_label(ax)
    plt.xlim(xlim)
    plt.subplot(2,1,2)    
    prh_hist(L0_protein_plot,normed=True,bins=bins,
             label="DNA+PRC2" + n_str(n_protein),**kw_protein)
    PlotUtilities.lazyLabel("L0 (microns)","P (1/microns)","")
    plt.xlim(xlim)
    PlotUtilities.savefig(fig,out_dir + "hist.png",
                          subplots_adjust=dict(hspace=0.03))
    for obj in objs_all:
        # plot each image with all the traces overlayed
        fig = PlotUtilities.figure()        
        plt.imshow(obj.image.height_nm_rel())
        # plot each DNA trace
        for o in obj.worm_objects:
            xy_abs = o.inf.x_y_abs
            color = "g" if o.has_dna_bound_protein else "r"
            plt.plot(*xy_abs,color=color,linewidth=0.25)
        out_name = os.path.basename(obj.image_path)
        PlotUtilities.savefig(fig,out_dir + out_name + ".png")



if __name__ == "__main__":
    run()
