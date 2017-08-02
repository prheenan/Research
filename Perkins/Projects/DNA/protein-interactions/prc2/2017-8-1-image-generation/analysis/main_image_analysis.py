# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,matplotlib,os,copy

sys.path.append("../../../../../../../../")

from GeneralUtil.python import GenUtilities,PlotUtilities

class worm_object:
    def __init__(self,x,y,text_file):
        self.x = x
        self.y = y
        self.L0_pixels = get_contour_length(self.x,self.y)
        self.L0_meters = None
        self.file_name = text_file
    def set_meters_per_pixel(self,m_per_pixel):
        self.m_per_pixel = m_per_pixel
        self.L0_meters = self.L0_pixels * self.m_per_pixel
    @property
    def has_dna_bound_protein(self):
        return "protein" in self.file_name.lower()
    @property
    def L0(self):
        assert self.L0_meters is not None
        return self.L0_meters       

class tagged_image:
    def __init__(self,image_path,worm_objects):
        self.image_path = image_path
        self.image = plt.imread(self.image_path)
        self.worm_objects = worm_objects
        cond =  [w.has_dna_bound_protein for w in self.worm_objects]
        self.protein_idx = np.where(cond)[0]
        self.dna_only_idx = np.where(~np.array(cond))[0]
    def subset(self,idx):
        return [copy.deepcopy(self.worm_objects[i]) for i in idx]
    def traces_dna_protein(self):
        return self.subset(self.protein_idx)
    def traces_dna_only(self):
        return self.subset(self.dna_only_idx)    
    def _L0(self,subset):
        return np.array([s.L0 for s in subset])
    def L0_protein_dna(self):
        return self._L0(self.traces_dna_protein())
    def L0_dna_only(self):
        return self._L0(self.traces_dna_only())
        

def get_x_and_y_arrays(text_file):
    data = np.loadtxt(text_file)
    return data[:,0],data[:,1]
    
def get_contour_length(x,y):
    return np.sum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
    
def get_x_y_and_contour_lengths(text_files):
    to_ret = []
    for t in text_files:
        x,y = get_x_and_y_arrays(t)
        to_ret.append(worm_object(x,y,t))
    return to_ret     
  
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "./in/"
    out_dir = "./out/"
    ext = ".tiff"
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
        these_text_files= [ t for t in text_files if str(file_no_ext) in str(t)]
        objs_tmp = get_x_y_and_contour_lengths(these_text_files)
        im = plt.imread(file_name)
        assert (im.shape[0] == im.shape[1]) 
        assert (im.shape[0] == size_images_pixels)
        # POST: dimensions are OK 
        objs_all.append(tagged_image(file_name,objs_tmp))
    # POST: we have all the data. Go ahead and set the conversion factor
    for o in objs_all:
        for wlc in o.worm_objects:
            wlc.set_meters_per_pixel(conversion_meters_per_px)
    # POST: all the contour lengths are set in 'real' units ]  
    L0_protein = np.concatenate([o.L0_protein_dna() for o in objs_all])
    L0_dna = np.concatenate([o.L0_dna_only() for o in objs_all])
    print(np.mean(L0_protein),np.std(L0_protein))
    print(np.mean(L0_dna),np.std(L0_dna))
    for obj in objs_all:
        fig = PlotUtilities.figure()        
        plt.imshow(obj.image)
        for o in obj.worm_objects:
            color = "g" if o.has_dna_bound_protein else "r"
            plt.plot(o.x,o.y,color=color,linewidth=0.25)
        out_name = os.path.basename(obj.image_path)
        PlotUtilities.savefig(fig,out_dir + out_name + ".png")




if __name__ == "__main__":
    run()
