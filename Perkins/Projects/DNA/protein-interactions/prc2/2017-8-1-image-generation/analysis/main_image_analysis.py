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
import re

from GeneralUtil.python import GenUtilities,PlotUtilities,CheckpointUtilities
from Research.Perkins.AnalysisUtil.Images import PolymerTracing,PolymerPlotting,\
    ImageUtil

def get_x_and_y_arrays(text_file):
    """
    Returns: the x and y columns (0 and 1) of text_file 
    """
    data = np.loadtxt(text_file)
    x = data[:,0]
    y = data[:,1]
    assert ((x >= 1) & (y > 1)).all()
    return x,y
    
def get_contour_length(x,y):
    """
    Returns: the contour length of the line L_i=(x,y)_i
    """
    return np.sum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
    
def get_x_y_and_contour_lengths(text_files):
    """
    Gets all the worm_objects associated with text_files
    
    Args:
        text_files: list of files, each has columns like <x>,<y>
    Returns:
        list of worm_object 
    """
    to_ret = []
    for t in text_files:
        x,y = get_x_and_y_arrays(t)
        to_ret.append(PolymerTracing.worm_object(x,y,t))
    return to_ret     

def prh_hist(data,**hist_kw):
    counts,_,_ = plt.hist(data,**hist_kw)
    y = max(counts * 1.05)
    ax = plt.gca()
    try:
        plt.boxplot(data,positions=[y],vert=False,manage_xticks=False,
                    meanline=True,showmeans=True,
                    flierprops=dict(color='k',markersize=1))
    except IndexError:
        # not enough data 
        pass
   
def print_info(dist,name):
    print("For {:s}:".format(name))
    # print in nm 
    # round to the tens place (pretty sure are errors are ~100nm)
    sanit = lambda x: int(np.round(x*1e9,-1))
    print("\t Median: {:.0f} nm".format(sanit(np.median(dist))))
    print("\t Mean  : {:.0f} nm".format(sanit(np.mean(dist))))
    print("\t Std   : {:.0f} nm".format(sanit(np.std(dist))))

def get_id(str_v):
    pat = r"""
           .+?    # anything (non greedy)
           Image  # the literal 'image'
           (\d{4})# 4 digits
           .?     # anything
           """
    match = re.match(pat,str_v,re.VERBOSE)
    assert match is not None
    id_v = match.group(1)
    return id_v
       
    
def yield_files(image_files,text_files,size_images_meters):
    text_image_ids = [get_id(t) for t in text_files ]
    for file_name in image_files:
        file_no_ext = file_name.replace(".pkl","")
        file_no_number = file_name.rsplit("_",1)[0]
        # get the text files that match 
        these_text_files= [text_files[i] for i,t in enumerate(text_image_ids)
                            if get_id(file_no_number) == str(t)]
        image_obj = CheckpointUtilities.lazy_load(file_name)
        # only look at images of size size_images_meters
        if ( abs(image_obj.range_meters -size_images_meters) > 1e-6):
            continue
        objs_tmp = get_x_y_and_contour_lengths(these_text_files)
        # POST: dimensions are OK 
        img = PolymerTracing.tagged_image(image_obj,objs_tmp,file_no_number)
        yield img
        
def read_images(in_dir,cache_dir):
    """
    reads in all images / pkl files / annotation traces from in_dir, saving
    intermediate results to cache_dir
    
    Args:
        in_dir: directory with pkl files and txt files as input to yield_files
        cache_dir: where to put the objects
    Returns:
        list of tagged_image objects 
    """
    ext = ".pkl"
    ext_text = ".txt"
    image_files = GenUtilities.getAllFiles(in_dir,ext=ext)
    text_files = GenUtilities.getAllFiles(in_dir,ext=ext_text)
    # filter only to those image files which are 2um
    size_images_meters = 2e-6
    size_images_pixels = 512
    conversion_meters_per_px = size_images_meters / size_images_pixels  
    objs_all = list(yield_files(image_files,text_files,size_images_meters))
    name_func = lambda i,d: "{:s}".format(d.file_name,i)
    load_func = lambda : yield_files(image_files,text_files,size_images_meters)
    objs_all = CheckpointUtilities.multi_load(cache_dir=cache_dir,
                                              load_func=load_func,
                                              force=False,name_func=name_func)
    for o in objs_all:
        im = o.image.height_nm_rel()
        assert (im.shape[0] == im.shape[1])
        assert (im.shape[0] == size_images_pixels)
    return objs_all

    
def make_contour_length_plot(objs_all,objs_0x,step_size_m=50e-9,ylim_max=None):
    L0_protein = np.concatenate([o.L0_protein_dna() for o in objs_all])
    L0_dna = np.concatenate([o.L0_dna_only() for o in objs_all])
    # only use the negative control if it exists
    # (note: we still plot it if so)
    has_negative_control = len(objs_0x) > 0
    if (has_negative_control):
        L0_dna_0x = np.concatenate([o.L0_dna_only() for o in objs_0x])
        # print off the stats. of we have them
        print_info(L0_dna_0x,"DNA+PRC2")        
    else:
        L0_dna_0x = np.array([])
    print_info(L0_dna,"DNA, not incubcated ")               
    print_info(L0_protein,"DNA")
    n_str = lambda n: "\n(N={:d})".format(n.size)
    sanit_L0 = lambda x: x*1e6
    L0_dna_plot = sanit_L0(L0_dna)
    L0_protein_plot = sanit_L0(L0_protein)
    L0_dna_0x = sanit_L0(L0_dna_0x)
    xmax = np.max(np.concatenate([L0_dna_plot,L0_protein_plot]))
    step = sanit_L0(step_size_m)
    bins = np.arange(start=0,stop=xmax+step,step=step)
    xlim = [0,xmax*1.1]
    legend_color = 'rebeccapurple'
    lazy_kwargs = dict(legend_kwargs=dict(color=legend_color),frameon=False,loc=(0.05,0.15))
    kw_dna = dict(color='g',alpha=0.3,lazy_kwargs=lazy_kwargs)
    kw_protein = dict(color='b',hatch='//',alpha=0.7,lazy_kwargs=lazy_kwargs)
    kw_dna_only = dict(color='r',alpha=0.3,hatch='x',lazy_kwargs=lazy_kwargs)
    ax0 = plt.subplot(3,1,1)
    color = 'rebeccapurple'
    plot_histogram(L0_dna_0x,bins,label="-/-" + n_str(L0_dna_0x),
                   **kw_dna_only)
    PlotUtilities.no_x_label(ax0)
    PlotUtilities.xlabel("")
    PlotUtilities.title("+/+ : DNA incubated with PRC2 / DNA bound to PRC2",
                        color=legend_color)
    ax1 = plt.subplot(3,1,2)
    plot_histogram(L0_dna_plot,bins,label="+/-" + n_str(L0_dna_plot),**kw_dna)

    ax2 = plt.subplot(3,1,3)
    plot_histogram(L0_protein_plot,bins,label="+/+" + n_str(L0_protein_plot),
                   **kw_protein)
    axs = [ax0,ax1,ax2]
    if (ylim_max is None):
        ylim_max_arr = [a.get_ylim()[1] for a in axs]
        ylim_max = max(ylim_max_arr) * 1.1
    ylim = [0,ylim_max]
    for ax_tmp in axs:
        ax_tmp.set_ylim(ylim)
        PlotUtilities.tom_ticks(ax=ax_tmp,num_major=1,change_x=False)
    for ax_tmp in axs[:-1]:
        PlotUtilities.no_x_label(ax_tmp)
        PlotUtilities.xlabel("",ax=ax_tmp)
    
def plot_histogram(data_plot,bins,lazy_kwargs=dict(),**kw_plot):
    """
    plots data_to_plot (y) histogrammed to bins. **kw_plot passed to prh_hist
    """
    micron_str = "$\mathrm{\mu m}$"
    prob_str = "$P$ (1/" + micron_str + ")"
    prh_hist(data_plot,normed=True,bins=bins,
             **kw_plot)
    PlotUtilities.lazyLabel("Contour length $L_0$ (" + micron_str + ")",prob_str,"",
                            **lazy_kwargs)
        
def plot_all_objects(out_dir,objs_all):
    """
    plots all the annotations on all the objects.
    """
    for obj in objs_all:
        # plot each image with all the traces overlayed
        fig = PlotUtilities.figure((3.5,4))
        ax1 = plt.subplot(2,1,1)
        plot_image(obj,ax=ax1,fig=fig,colorbar_kw=dict(add_space_only=True))
        PlotUtilities.no_x_label(ax1)
        ax2 = plt.subplot(2,1,2)
        plot_annotated_object(obj,ax=ax2,fig=fig)
        out_name = os.path.basename(obj.image_path)        
        PlotUtilities.savefig(fig,out_dir + out_name + ".png",
                              subplots_adjust=dict(hspace=0.05))    

def plot_image(obj,ax,fig,colorbar_kw=dict()):
    """
    Plots a single image, relative to the median (assumed to be the surface)
    
    Args:
        obj: tagged_image object
        ax: the axis to plot on
        fig: the figur to apply upon
        colorbar_kw: passed to ImageUtil.smart_colorbar
    Returns:
        nothing
    """
    imshow_kw = dict(vmin=0,vmax=1.25)                          
    height_nm_rel_surface = obj.image.height_nm_rel()
    height_nm_rel_surface -= np.median(height_nm_rel_surface)
    im = plt.imshow(height_nm_rel_surface,cmap=plt.cm.Greys_r,**imshow_kw)
    ImageUtil.smart_colorbar(im,ax=ax,fig=fig,**colorbar_kw)
    
def plot_annotated_object(obj,ax,fig):
    """
    plots an object with its annotations. See: plot_image
    """
    plot_image(obj,ax,fig)
    # plot each DNA trace
    for o in obj.worm_objects:
        xy_abs = o.inf.x_y_abs
        color = "g" if o.has_dna_bound_protein else "r"
        ax.plot(*xy_abs,color=color,linewidth=0.2)

def detailed_plot(in_dir,in_dir_0x,out_dir,cache_base):
    GenUtilities.ensureDirExists(out_dir)    
    objs_0x = read_images(in_dir_0x,cache_dir=cache_base + "0x/")    
    objs_all = read_images(in_dir,cache_dir=cache_base + "1x/")
    print(objs_all)
    kw = dict(min_m = 0,max_m = 125e-9)
    polymer_info_obj = PolymerTracing.ensemble_polymer_info(objs_all,**kw)
    fig = PlotUtilities.figure()
    PolymerPlotting.plot_angle_information(polymer_info_obj)
    PlotUtilities.savefig(fig,out_dir + "angles.png")
    # POST: all the contour lengths are set in 'real' units ]
    fig = PlotUtilities.figure()
    make_contour_length_plot(objs_all,objs_0x)
    PlotUtilities.savefig(fig,out_dir + "2017-10-4-histograms.png",
                          subplots_adjust=dict(hspace=0.07))    

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "in/"
    mM_fmt = lambda mM : "{:d}mM".format(mM)
    args_mM = lambda mM : \
        dict(in_dir="./in/in_{:s}/".format(mM_fmt(mM)),
             out_dir="./out/out_{:s}/".format(mM_fmt(mM)),
             cache_base="./cache/cache_{:s}".format(mM_fmt(mM)))
    detailed_plot(in_dir_0x="./in/in_empty/",**args_mM(120))  
    detailed_plot(in_dir_0x="./in/in_empty/",**args_mM(10))
    #plot_all_objects(out_dir,objs_all)        


if __name__ == "__main__":
    run()
