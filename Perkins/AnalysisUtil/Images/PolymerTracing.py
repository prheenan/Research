# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,scipy,copy

from scipy.interpolate import splprep, splev, interp1d,UnivariateSpline

from scipy.stats import binned_statistic
        
from skimage.morphology import skeletonize,medial_axis,dilation
from skimage import measure
from scipy.interpolate import LSQUnivariateSpline
from skimage.segmentation import active_contour
from skimage.filters import gaussian


class spline_info(object):
    # holds low-level information on the actual fitting used. Everything
    # is in units of pixels
    def __init__(self, u, tck, spline, deriv,x0_px,y0_px):
        self.u = u
        self.tck = tck
        self.spline = spline
        self.deriv = deriv
        self.x0_px = x0_px
        self.y0_px = y0_px

class polymer_info(object):
    # holds low-level information about the polymer itself
    def __init__(self,L_m,cos_angle,Lp_m,L0_m):
        self.L_m = L_m
        self.cos_angle = cos_angle
        self.Lp_m = Lp_m
        self.L0_m = L0_m

class spline_fit_obj(object):
    # spline_fit_obj: class for holding all information about a given spline fit
    def __init__(self,image_cropped,image_threshold,polymer_info_obj,
                 fit_spline):
        self.image_cropped=image_cropped,
        self.image_threshold=image_threshold
        self.polymer_info_obj = polymer_info_obj
        self.fit_spline= fit_spline
    @property
    def x0_y0(self):
        return self.fit_spline.x0_px,self.fit_spline.y0_px
    @property
    def x_y_rel(self):
        return self.fit_spline.spline
    @property
    def x_y_abs(self):
        """
        Returns: the x and y coordinates in the absolute, original image coords
        """
        x_rel,y_rel = self.x_y_rel
        x0,y0 = self.x0_y0
        return x_rel+x0,y_rel+y0
    @property
    def L0_m(self):
        return self.polymer_info_obj.L0_m
    @property
    def Lp_m(self):
        return self.polymer_info_obj.Lp_m

class worm_object(object):
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
        self.file_name = text_file
        self.inf = None
    def assert_fit(self):
        assert self.inf is not None
    def set_spline_info(self,inf):
        self.inf = inf
    @property
    def has_dna_bound_protein(self):
        """
        Returns: True iff this trace was tagged as having protein 
        """
        return "protein" in self.file_name.lower()
    @property
    def Lp(self):
        """
        Returns: the persistence length, in meters
        """
        self.assert_fit()    
        return (self.inf.Lp_m)
    @property
    def L0(self):
        """
        Returns: the contour length, in meters
        """
        self.assert_fit()
        return (self.inf.L0_m)

class tagged_image:
    def __init__(self,image,worm_objects,file_no_number):
        """
        Grouping of an images and the associated traces on items
        
        Args:
            image_path: the file name of the image
            worm_objects: list of worm_objects associated with this image
        """
        self.image_path = file_no_number
        self.image = image
        self.worm_objects = worm_objects
        cond =  [w.has_dna_bound_protein for w in self.worm_objects]
        assert len(cond) > 0 , "{:s} has no tagged data".format(image.Meta.Name)
        for w in worm_objects:
            tmp_fit = spline_fit(self.image,x=w._x_raw,y=w._y_raw)
            w.set_spline_info(tmp_fit)
        # POST: as least something to look at 
        self.protein_idx = np.where(cond)[0]
        self.dna_only_idx = np.where(~np.array(cond))[0]
    @property
    def file_name(self):
        return self.image_path.rsplit("/",1)[-1]
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
    def _subset_f(self,f,subset):
        return np.array([f(s) for s in subset])
    def _L0(self,subset):
        """
        Returns the contour length of the given subset 
        """
        return self._subset_f(lambda s: s.L0,subset)
    def _Lp(self,subset):
        """
        Returns the persistence length of the given subset 
        """
        return self._subset_f(lambda s: s.Lp,subset)
    def _f_dna_protein(self,f):
        return f(self.traces_dna_protein())
    def _f_dna_only(self,f):
        return f(self.traces_dna_only())
    def L0_protein_dna(self):
        """
        Returns: the contour length, in meters, of the traces on DNA-protein 
        complex in this image 
        """    
        return self._f_dna_protein(self._L0)
    def L0_dna_only(self):
        """
        Returns: the contour length, in meters, of the traces on DNA in this
        image 
        """
        return self._f_dna_only(self._L0)
        

def crop_slice(data,f=1):
    v_min,v_max = min(data),max(data)
    v_range = v_max-v_min
    n_v = int(np.ceil(f * v_range))
    v_low,v_high = max(0,v_min-n_v),v_max+n_v
    return slice(int(v_low),int(v_high),1)

def get_region_of_interest(height_cropped_nm,background_image,threshold_nm=0.0):
    """
    Returns: single-region of interest
    
    Args:
        height_cropped_nm: image, elements are in nm
        background_image: value considered the background for height_cropped_nm
        thresold_nm: for thresholding, the minimum above the background
        to be considered not noise

    Returns:
        image, same shape as height_cropped_nm. everything not in the largest
        skeleton region is set to zero
    """
    # threshold anything less than x nM
    image_thresh = height_cropped_nm.copy()
    image_thresh[np.where(image_thresh < background_image+threshold_nm)] = 0 
    # binarize the image
    image_binary = image_thresh.copy()
    image_binary[np.where(image_binary > 0)] = 1
    # skeletonize the image
    image_skeleton = skeletonize(image_binary)
    image_label = measure.label(image_skeleton,background=0)
    props = measure.regionprops(image_label) 
    diameters = [p.equivalent_diameter for p in props]
    # XXX just use x,y, dialyze that?
    assert len(diameters) > 0 , "Couldn't find any objects in region"
    max_prop_idx = np.argmax(diameters)
    largest_skeleton_props = props[max_prop_idx]
    # zero out everything not the one we want 
    skeleton_zeroed = np.zeros(image_thresh.shape)
    # take the largest object in the view, zero everything else
    # order the points in that object by the x-y point 
    x_skel = largest_skeleton_props.coords[:,0]
    y_skel = largest_skeleton_props.coords[:,1]
    for x_tmp,y_tmp in zip(x_skel,y_skel):
        skeleton_zeroed[x_tmp,y_tmp] = 1
    # dilated skeleton; make it 'fat'
    dilation_size = 3
    selem = np.ones((dilation_size,dilation_size))
    skeleton_zeroed = dilation(skeleton_zeroed,selem=selem)
    # mask the original data with the skeletonized one
    image_single_region = skeleton_zeroed * height_cropped_nm
    return image_single_region

def get_L_and_mean_angle(cos_angle,L,n_bins,min_cos_angle = np.exp(-2)):
    """
    Gets L and <cos(theta(L))>

    Args:
        cos_angle: length N, element i is angle between two segmens
        L: length N, element i is contour length between same segments as 
        cos_angle
        
        n_bins: we will average cos_angle in this many bins from its min to max
       
        min_cos_angle: we cant use when <cos(Theta)> <= 0,since that would go
        negative when we take a log. So, only look where above this value

    Returns:
       tuple of L_[avg,j],<Cos(Theta_[avg,j])>, where j runs 0 to n_bins-1
    """
    bins = np.linspace(min(cos_angle),max(L),num=n_bins)
    mean_cos_angle,edges,_ = binned_statistic(x=L,values=cos_angle,
                                              bins=bins)
    # last edge is right bin
    edges = edges[:-1]
    # filter to the bins with at least f% of the total size
    f_min_size = 1/(bins.size)
    values,_ = np.histogram(a=L,bins=bins)
    bins_with_data = np.where(values > 0)
    mean_cos_angle = mean_cos_angle[bins_with_data]
    edges = edges[bins_with_data]
    # only look at where cos(theta) is reasonable positive, otherwise we 
    # cant take a log. This amounts to only looking in the upper quad
    good_idx = np.where((mean_cos_angle > min_cos_angle))
    sanit = lambda x: x[good_idx]
    mean_cos_angle = sanit(mean_cos_angle)
    edges = sanit(edges)
    assert edges.size > 0 , "Couldn't find data to fit"
    return edges,mean_cos_angle

def Lp_log_mean_angle_and_coeffs(L,mean_cos_angle):
    """
    Returns: the persistence length, -Log(<Cos(Theta(L))), and linear polynomial
             coefficients for a given <Cos(Theta(L))>
    """
    log_mean_angle = -np.log(mean_cos_angle)
    # fit to -log<Cos(angle)> to edges_nm
    coeffs = np.polyfit(x=L,y=log_mean_angle,deg=1)
    persistence_length = 1/coeffs[0]
    return persistence_length,log_mean_angle,coeffs

def spline_fit(image_obj,x,y):
    image = image_obj.height_nm_rel()
    slice_x = crop_slice(x)
    slice_y = crop_slice(y)
    # the matrix is transposed, so swap x and y
    background_image = np.percentile(image,50)
    image_cropped = image[slice_y,slice_x]
    x0,y0 = slice_x.start,slice_y.start
    x_rel = x-x0
    y_rel = y-y0
    image_single_region = \
            get_region_of_interest(height_cropped_nm=image_cropped,
                                   background_image=background_image)
    good_idx = np.where(image_single_region > 0)
    xy_rel = np.array((x_rel,y_rel)).T
    # fit a spline to r(t)=(x(t),y(t)) to the manually-tagged data
    u,tck,spline,deriv = _u_tck_spline_and_derivative(x_rel,y_rel)
    assert image.shape[0] == image.shape[1] , \
        "Non square image unsupported"
    m_per_px = (image_obj.range_meters/image.shape[0])
    cos_angle,flat_L0,L0_px = \
        angles_and_contour_lengths(spline,deriv,
                                   min_change_px=0,
                                   max_change_px=100e-9/m_per_px)
    # POST: cos_angle and flat_L0 and reasonable
    edges,mean_cos_angle =  get_L_and_mean_angle(cos_angle,flat_L0,n_bins=50)
    L_m = edges * m_per_px
    L0_m = L0_px * m_per_px
    Lp_m,log_mean_angle,_ = \
        Lp_log_mean_angle_and_coeffs(L_m,mean_cos_angle)
    # do some data checking.        
    assert L0_m > 0 , "L0 must be positive"
    assert Lp_m > 0 , "Lp must be positive"
    # POST: most basic polymer stuff is OK.
    fit_spline_info = spline_info(u,tck,spline,deriv,x0_px=x0, y0_px=y0)
    polymer_info_obj =polymer_info(L_m=L_m, cos_angle=cos_angle,
                                   Lp_m=Lp_m,L0_m=L0_m)
    return  spline_fit_obj(image_cropped=image_cropped,
                           image_threshold=image_single_region,
                           polymer_info_obj=polymer_info_obj,
                           fit_spline=fit_spline_info)


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
    L0 = contour_lengths[-1]
    n = x_spline.shape[0]
    contour_length_matrix = _difference_matrix(contour_lengths,contour_lengths)
    cos_angle_matrix = _dot_matrix(deriv_unit_vector.T,deriv_unit_vector.T)
    cos_angle_matrix = np.maximum(-1,cos_angle_matrix)
    cos_angle_matrix = np.minimum(1,cos_angle_matrix)
    # only look at the upper triangular part
    idx_upper_tri = np.triu_indices(n)
    idx_upper_tri_no_diag =np.triu_indices(n,k=1)
    # angle should be finite, between -1 and 1, also be finite
    assert ((cos_angle_matrix >= -1) & \
            (cos_angle_matrix <= 1) & \
            (np.isfinite(cos_angle_matrix))).all()
    # upper diagonal should have >0 contour length
    assert (contour_length_matrix[idx_upper_tri_no_diag] > 0).all() , \
        "Contour lengths should be positive"
    # POST: contour lengths and angles make sense; we only want upper triangular
    # (*including* the trivial 0,0 point along the diagonal)
    contour_length_matrix = contour_length_matrix[idx_upper_tri]
    cos_angle_matrix = cos_angle_matrix[idx_upper_tri]
    # POST: matrix is filled in, determine where the value are valid
    ok_idx = np.where( (contour_length_matrix > min_change_px) &
                       (contour_length_matrix < max_change_px))
    sanit = lambda x: x[ok_idx].flatten()
    flat_L0 = sanit(contour_length_matrix)
    sort_idx = np.argsort(flat_L0)
    # sort everything by L0...
    flat_L0 = flat_L0[sort_idx]
    cos_angle = sanit(cos_angle_matrix)[sort_idx]
    return cos_angle,flat_L0,L0

def _spline_u_and_tck(x,y,k=3,s=None,num=None):
    """
    fits a line r(u)=[x(u),y(u)], where u is determined implicitly

    Args:
        x,y: the coordinates to fit to
        k: the degree of the spline
        s: the smoothign factor, default (0) is no smoothing
    Returns:
        tuple of <u,tck>, where tck is needed for (e.g.) splev
    """
    if (num is None):
        num = len(x) * 10
    if (s is None):
        s = 0
    tck,_ = splprep([x,y],per=0,
                    s=s,k=k,u=None)
    u = np.linspace(0,1,num=num,endpoint=True)
    return u,tck

def _u_tck_spline_and_derivative(x,y,*kw):
    """
    Args:
        see _spline_u_and_tck
    Returns:
        see _spline_u_and_tck, plus the evaluated spline and derivative
        (as an in-order tuple)
    """
    u,tck = _spline_u_and_tck(x,y)
    spline = splev(x=u,tck=tck,der=0)
    deriv = splev(x=u,tck=tck,der=1)
    return u,tck,spline,deriv

def _difference_matrix(v1,v2):
    """
    Args:
        v<1/2>: the two vectors to subtract, size n and m
    Returns:
        matrix, size n X m, element i,j  is v1[j]-v2[i]
    """
    return (v2 - v1[:,np.newaxis])

def _dot_matrix(v1,v2):
    """
    Args:
        see _difference_matrix
    Returns:
        matrix M, where element i,j is v1[i] . v2[j]
    """
    return np.dot(v1,v2.T)
