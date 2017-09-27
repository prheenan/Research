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
from scipy.interpolate import splprep, splev, interp1d,UnivariateSpline

from scipy.stats import binned_statistic

class worm_object:
    def __init__(self,x,y,text_file,spline_kwargs=dict(k=3)):
        """
        object for keeping track of an x,y trace
        
        Args:
            x,y: thew coordinates 
            text_file: where this trace came from 
        """
        self._x_raw = x
        self._y_raw = y
        self.spline_kwargs=spline_kwargs
        # XXX fix...
        #self.x,self.y = spline_fit(x,y,**spline_kwargs)
        #self.L0_pixels = get_contour_length(self.x,self.y)
        self.L0_meters = None
        self.file_name = text_file
    def set_meters_per_pixel(self,m_per_pixel):
        """
        Sets the scale in meters/pixel to m_per_pixel. Needed for L0_meters
        """
        self.m_per_pixel = m_per_pixel
        self.L0_meters = self.L0_pixels * self.m_per_pixel
    @property
    def has_dna_bound_protein(self):
        """
        Returns: True iff this trace was tagged as having protein 
        """
        return "protein" in self.file_name.lower()
    @property
    def L0(self):
        assert self.L0_meters is not None
        return self.L0_meters       

class tagged_image:
    def __init__(self,image,worm_objects):
        """
        Grouping of an images and the associated traces on items
        
        Args:
            image_path: the file name of the image
            worm_objects: list of worm_objects associated with this image
        """
        self.image = image
        self.worm_objects = worm_objects
        cond =  [w.has_dna_bound_protein for w in self.worm_objects]
        assert len(cond) > 0 , "{:s} has no tagged data".format(image.Meta.Name)
        for w in worm_objects:
            spline_fit(self.image,w)
        # POST: as least something to look at 
        self.protein_idx = np.where(cond)[0]
        self.dna_only_idx = np.where(~np.array(cond))[0]
    def subset(self,idx):
        return [copy.deepcopy(self.worm_objects[i]) for i in idx]
    def traces_dna_protein(self):
        """
        Returns: the traces on this image of protein bound to DNA 
        """
        return self.subset(self.protein_idx)
    def traces_dna_only(self):
        """
        Returns: the traces on this object on DNA (only)
        """
        return self.subset(self.dna_only_idx)    
    def _L0(self,subset):
        """
        Returns the contour length of the given subset 
        """
        return np.array([s.L0 for s in subset])
    def L0_protein_dna(self):
        """
        Returns: the contour length, in meters, of the traces on DNA-protein 
        complex in this image 
        """    
        return self._L0(self.traces_dna_protein())
    def L0_dna_only(self):
        """
        Returns: the contour length, in meters, of the traces on DNA in this
        image 
        """
        return self._L0(self.traces_dna_only())
        
def crop_slice(data,f=0.3):
    v_min,v_max = min(data),max(data)
    v_range = v_max-v_min
    n_v = int(np.ceil(f * v_range))
    v_low,v_high = max(0,v_min-n_v),v_max+n_v
    return slice(int(v_low),int(v_high),1)
        
from skimage.morphology import skeletonize,medial_axis,dilation
from skimage import measure
from scipy.interpolate import LSQUnivariateSpline
from skimage.segmentation import active_contour
from skimage.filters import gaussian

def spline_fit(image_obj,worm_object):
    image = image_obj.height_nm_rel()
    x = worm_object._x_raw
    y = worm_object._y_raw
    slice_x = crop_slice(x)
    slice_y = crop_slice(y)
    # the matrix is transposed, so swap x and y
    background_image = np.percentile(image,50)
    image_cropped = image[slice_y,slice_x]
    x_rel = x-slice_x.start
    y_rel = y-slice_y.start
    # threshold anything less than 0.2nM
    image_thresh = image_cropped.copy()
    image_thresh[np.where(image_thresh < background_image+0.2)] = 0 
    # binarize the image
    image_binary = image_thresh.copy()
    image_binary[np.where(image_binary > 0)] = 1
    # skeletonize the image
    image_skeleton = skeletonize(image_binary)
    image_label = measure.label(image_skeleton,background=0)
    props = measure.regionprops(image_label) 
    diameters = [p.equivalent_diameter for p in props]
    max_prop_idx = np.argmax(diameters)
    largest_skeleton_props = props[max_prop_idx]
    # zero out everything not the one we want 
    skeleton_zeroed = np.zeros(image_cropped.shape)
    # take the largest object in the view, zero everything else
    # order the points in that object by the x-y point 
    x_skel = largest_skeleton_props.coords[:,1]
    y_skel = largest_skeleton_props.coords[:,0]
    for x_tmp,y_tmp in zip(x_skel,y_skel):
        skeleton_zeroed[x_tmp,y_tmp] = 1
    # dilated skeleton; make it 'fat'
    dilation_size = 3
    selem = np.ones((dilation_size,dilation_size))
    skeleton_zeroed = dilation(skeleton_zeroed,selem=selem)
    # mask the original data with the skeletonized one
    image_single_region = skeleton_zeroed * image_cropped
    good_idx = np.where(image_single_region > 0)
    x_region,y_region = good_idx
    z_region = image_single_region[good_idx]
    xy_skel = np.array((x_region,y_region)).T
    xy_rel = np.array((x_rel,y_rel)).T
    # get the closest trace index (ie: 'hand picked' pixel')
    # to each pixel in the data
    distance_matrix = scipy.spatial.distance_matrix(xy_skel,xy_rel)
    closest_idx = np.argmin(distance_matrix,axis=1)
    # fit a spline to r(t)=(x(t),y(t)) to the manually-tagged data
    spline,deriv = _spline_and_derivative(x_rel,y_rel)
    x_spline,y_spline = spline
    nm_per_px = (2e-6/512) * 1e9
    cosine_of_angle,flat_L0 = \
        angles_and_contour_lengths(spline,deriv,
                                   min_change_px=0,max_change_px=100/nm_per_px)
    # do some checks to make sure the data are sensible
    assert ((cosine_of_angle <= 1) & (cosine_of_angle >= -1)).all()
    # POST: data are reasonable
    bins = np.linspace(0,max(flat_L0),num=50)
    mean_cos_angle,edges,_ = \
        binned_statistic(x=flat_L0,values=cosine_of_angle,bins=bins)
    # last edge is right bin
    edges = edges[:-1]
    # filter to the bins with at least f% of the total size
    f_min_size = 1/(1.5*bins.size)
    values,_ = np.histogram(a=flat_L0,bins=bins)
    # only look at where cos(theta) is reasonable positive, otherwise we 
    # cant take a log. This amounts to only looking in the upper quarant 
    good_idx = np.where(mean_cos_angle > np.exp(-2))
    sanit = lambda x: x[good_idx]
    mean_cos_angle = sanit(mean_cos_angle)
    edges = sanit(edges)
    edges_nm = edges*nm_per_px
    log_mean_angle = -np.log(mean_cos_angle)
    # fit to -log<Cos(angle)> to edges_nm
    coeffs = np.polyfit(x=edges_nm,y=log_mean_angle,deg=1)
    persistence_length_nm,offset_nm = 1/coeffs,coeffs[1]
    predicted = np.polyval(coeffs,x=edges_nm)
    plt.plot(edges_nm,log_mean_angle,'go',linewidth=0.5)
    plt.plot(edges_nm,predicted,'b--',linewidth=2)
    plt.show()

def angles_and_contour_lengths(spline,deriv,
                               min_change_px=0,max_change_px=np.inf):
    """
    gets Cos(Theta(i)) and L(i), where i runs along the spline order given,
    and L is the contour length between segments chosen at index i

    Args:
        spline: tuple of x_spline,y_spline -- x and y values of the line, size N
        deriv: the continuous derivative of spline, size N
        <min/max>_change_px: the minimum and maximum pixel changes
    Returns:
        tuple of Cos(Theta(i)) and L(i) 
    """
    # get the x and y coordinates of the spline
    x_spline,y_spline = spline
    x_deriv,y_deriv = deriv
    deriv_unit_vector = np.array((x_deriv,y_deriv))
    deriv_unit_vector /= np.sqrt(np.sum(np.abs(deriv_unit_vector**2),axis=0))
    assert ((np.sum(deriv_unit_vector**2,axis=0) -1) < 1e-6).all() , \
        "Unit vectors not correct"
    # POST: unit vector are normalized, |v| = 1
    dx_spline = np.array([0] + list(np.diff(x_spline)))
    dy_spline = np.array([0] + list(np.diff(y_spline)))
    # d_spline(i) is the change from i-i to i (zero if i=0)
    d_spline = np.sqrt(dx_spline**2 + dy_spline**2)
    assert (dx_spline <= d_spline).all()
    contour_lengths = np.cumsum(d_spline)
    n = x_spline.shape[0]
    contour_length_matrix = np.zeros((n,n))
    cos_angle_matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            # only fill in the upper part of the matrix; avoid double-counting
            contour_length_matrix[i,j] = contour_lengths[j] - contour_lengths[i]
            # cos(theta_[a,b]) = (a . b)/(|a|*|b|) = (a . b)
            # (for |a|=|b|=1)
            cos_angle_tmp = np.dot(deriv_unit_vector[:,i],
                                   deriv_unit_vector[:,j])
            cos_angle_matrix[i,j] = cos_angle_tmp
    # POST: matrix is fileld in, determine where the value are valid
    ok_idx = np.where( (contour_length_matrix > min_change_px) &
                       (contour_length_matrix < max_change_px))
    sanit = lambda x: x[ok_idx].flatten()
    flat_L0 = sanit(contour_length_matrix)
    sort_idx = np.argsort(flat_L0)
    # sort everything by L0...
    flat_L0 = flat_L0[sort_idx]
    cos_angle = sanit(cos_angle_matrix)[sort_idx]
    return cos_angle,flat_L0

def _spline_u_and_tck(x,y,k=3,s=None,num=None):
    if (num is None):
        num = len(x) * 20
    if (s is None):
        s = 0
    tck,_ = splprep([x,y],per=0,
                    s=s,k=k,u=None)
    u = np.linspace(0,1,num=num,endpoint=True)
    return u,tck

def _spline_and_derivative(x,y,*kw):
    u,tck = _spline_u_and_tck(x,y)
    spline = splev(x=u,tck=tck,der=0)
    deriv = splev(x=u,tck=tck,der=1)
    return spline,deriv

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
        to_ret.append(worm_object(x,y,t))
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
        objs_all.append(tagged_image(image_obj,objs_tmp))
    # POST: we have all the data. Go ahead and set the conversion factor
    for o in objs_all:
        for wlc in o.worm_objects:
            wlc.set_meters_per_pixel(conversion_meters_per_px)
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
        plt.imshow(obj.image)
        # plot each DNA trace
        for o in obj.worm_objects:
            color = "g" if o.has_dna_bound_protein else "r"
            plt.plot(o.x,o.y,color=color,linewidth=0.25)
        out_name = os.path.basename(obj.image_path)
        PlotUtilities.savefig(fig,out_dir + out_name + ".png")



if __name__ == "__main__":
    run()
