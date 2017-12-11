# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,re
sys.path.append("../../../../../../../../../../")
from GeneralUtil.python import GenUtilities,CheckpointUtilities,PlotUtilities
from Research.Perkins.AnalysisUtil.Images import PolymerTracing,\
    PolymerPlotting,ImageUtil
import argparse

def get_directory_command_line(description=""):
    parser = argparse.ArgumentParser(description)
    parser.add_argument("--input",type=str,metavar="input",
                        help="Please list the input directory")
    return parser.parse_args().input

def _dir_sanit(d):
    return d.replace("//","/")
    
def data_dir(base):
    return _dir_sanit("{:s}/pkl/".format(base))
    
def cache_dir(base,n,name):
    return _dir_sanit("{:s}/cache_{:d}_{:s}/".format(base,n,name))

def _traces_dir(base):
    return cache_dir(base,0,"spline_traces")
    
def _ensemble_dir(base):
    return cache_dir(base,1,"ensemble")
    
def _plot_dir(base):
    to_ret = _dir_sanit("{:s}/plot/".format(base))
    GenUtilities.ensureDirExists(to_ret)
    return to_ret
        
        
class EnsembleObject:
    def __init__(self,dna_only,dna_plus_protein,multimer=None,unknown=None):
        self.dna_only = dna_only
        self.dna_plus_protein = dna_plus_protein
        self.multimer = multimer
        self.unknown = unknown
        
def get_x_and_y_arrays(text_file):
    """
    Returns: the x and y columns (0 and 1) of text_file 
    """
    data = np.loadtxt(text_file,delimiter=" , ",comments="#")
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
        # read the first line (for possible protein locations)
        locations = []
        with open(t) as f:
            first_line = f.readline()
            idx_x_y_pattern = """<(\d+),(\d+),(\d+)>"""
            matches = re.finditer(idx_x_y_pattern,first_line, re.VERBOSE)
            locations = [int(i) for m in matches for i in m.groups()]
        worm_objs = PolymerTracing.worm_object(x,y,t)
        worm_objs.protein_locations = locations
        to_ret.append(worm_objs)
    return to_ret     

def prh_hist(data,y=None,**hist_kw):
    counts,_,_ = plt.hist(data,**hist_kw)
    if (y is None):
        y = max(counts * 1.2)
    ax = plt.gca()
    try:
        plt.boxplot(data,positions=[y],vert=False,manage_xticks=False,
                    meanline=True,showmeans=True,widths=y*0.07,
                    whis=[5,95],
                    boxprops=dict(linewidth=0.5),
                    meanprops=dict(linestyle='--'),
                    whiskerprops=dict(linewidth=0.5),
                    flierprops=dict(color='k',markersize=0.5,
                                    linewidth=0.25))
    except IndexError:
        # not enough data 
        pass
   
def print_info(dist,name):
    print("For {:s}:".format(name))
    # print in nm 
    # round to the tens place (pretty sure are errors are ~100nm)
    sanit = lambda x: int(np.round(x*1e9,-1)) if np.isfinite(x) else np.nan
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
        # get the text files that match
        m_id = get_id(file_name)
        these_text_files= [text_files[i] for i,t in enumerate(text_image_ids)
                            if m_id == str(t)]
        image_obj = CheckpointUtilities.lazy_load(file_name)
        # only look at images of size size_images_meters
        if ( abs(image_obj.range_meters -size_images_meters) > 1e-6):
            continue
        objs_tmp = get_x_y_and_contour_lengths(these_text_files)
        # POST: dimensions are OK 
        img = PolymerTracing.tagged_image(image_obj,objs_tmp,file_name)
        if (len(img.worm_objects) == 0):
            print("Warning, {:s} has no tagged data...".format(file_name))
            continue
        yield img
        
def read_images(in_dir,cache_dir,force=False):
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
    name_func = lambda i,d: "{:s}".format(d.file_name,i)
    load_func = lambda : yield_files(image_files,text_files,size_images_meters)
    objs_all = CheckpointUtilities.multi_load(cache_dir=cache_dir,
                                              load_func=load_func,
                                              force=force,name_func=name_func)
    for o in objs_all:
        im = o.image.height_nm_rel()
        assert (im.shape[0] == im.shape[1])
        assert (im.shape[0] == size_images_pixels)
    return objs_all        
