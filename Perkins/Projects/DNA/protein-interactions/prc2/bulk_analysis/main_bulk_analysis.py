# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,copy,scipy
sys.path.append("../../../../../../../")
from Research.Perkins.AnalysisUtil.Images import  ImageUtil
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from GeneralUtil.python import GenUtilities,PlotUtilities,CheckpointUtilities
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.filters import gaussian
from skimage.morphology import skeletonize
from skimage import measure,img_as_uint
from skimage.segmentation import active_contour
import networkx as nx

class transform:
    def __init__(self,name,function,imshow_kw=dict(cmap=plt.cm.afmhot)):
        self.name = name
        self.function = function
        self.imshow_kw = imshow_kw

def subtract_background(image,deg=2,**kwargs):
    """
    subtracts a line of <deg> from each row in <images>
    """
    image = image.height
    to_ret = image.copy().T
    shape = image.shape
    coords = np.arange(shape[1])
    coeffs = np.array(np.polyfit(x=coords,y=image,deg=deg,**kwargs))
    n_rows = shape[0]
    for i in range(n_rows):
        to_ret[i,:] -= np.polyval(coeffs[:,i],x=coords)
    return to_ret.T

def realistic_min_max(image,q_low=90,q_high=99.9):
    data = image.height_nm()
    low,high = np.percentile(data.flatten(),[q_low,q_high])
    return low,high
    
def plot_with_background_corrected(args,imshow_kw_list=None):
    n = len(args)
    fig = PlotUtilities.figure((3.5,2*n))
    for i,a in enumerate(args):
        ax = plt.subplot(n,1,(i+1))
        m_list = imshow_kw_list[i]
        if (m_list['cmap'] == plt.cm.afmhot):
            vmin,vmax = realistic_min_max(a)
        else:
            vmin,vmax = 0,None
        imshow_kwargs = dict(vmin=vmin,vmax=vmax,**m_list)
        im = ImageUtil.make_image_plot(a,pct=50,imshow_kwargs=imshow_kwargs)
        if (i == 0):
            ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig)
        if (i < n-1):
            PlotUtilities.xlabel("",ax=ax)
            PlotUtilities.no_x_label(ax=ax)
            ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig,add_space_only=True)
        else:
            ImageUtil.smart_colorbar(im=im,ax=ax,fig=fig,add_space_only=True)
    return fig


def _safe_apply(images,f):
    to_ret = []
    for ex in images:
        tmp = copy.deepcopy(ex)
        ret = f(ex)
        assert type(ret) is np.ndarray , "{:s} didn't return array".format(f)
        tmp.height = ret
        to_ret.append(tmp) 
    return to_ret 

def threshold(im,threshold_nm,rel_pct=50):
    """
    Returns: im, height zeroed where height-rel_pct_of_height < threshold_nm
    """
    height_rel = im.height_nm() 
    height_rel -= np.percentile(height_rel,rel_pct)
    zero_idx = np.where(height_rel < threshold_nm)
    height_new = copy.deepcopy(im.height)
    height_new[zero_idx] = 0
    return height_new

def binarize(image):
    """
    binarizes a single image: set to 1 where the image is non-zero
    """
    binary = copy.deepcopy(image)
    binary[binary > 0] = 1
    return binary

def correct_background(images,**kw):
    """
    See: threshold_images, except subtracts the AFM background 
    """
    return _safe_apply(images,lambda x: subtract_background(x,**kw))


def blur_images(images,sigma=0.5,**kw):
    """
    See: thresholdimages, except adds a gaussian blur with sigma
    """
    return _safe_apply(images,lambda x: gaussian(x.height,sigma=sigma,**kw))

def threshold_images(images,threshold_nm=0.2):
    """
    thresholds the heights for each of images. pass by copy
    """
    return _safe_apply(images,lambda x: threshold(x,threshold_nm))

def binarize_images(images):
    """
    binarizes the heights for each of images. pass by copy
    """
    return _safe_apply(images,lambda x: binarize(x.height))

def skeletonize_images(images):
    """
    returns: the skeletonized version of the (assumed already binary) images
    """
    return _safe_apply(images,lambda x: skeletonize(x.height))

def label_images(images):
    """
    returns: the labelled versions (ie: connected components) of the (assumed
             skeletonized) images
    """
    return _safe_apply(images,lambda x: measure.label(x.height,background=0))

def skeleton_filter(images):
    to_ret = []
    for i in images:
        tmp = copy.deepcopy(i)
        props = measure.regionprops(tmp.height)
        n = len(props)
        n_lost = 0
        for p in props:
            diameter = p.equivalent_diameter
            if (diameter < 5):
                for i,j in p.coords:
                    tmp.height[i,j] = 0 
                n_lost += 1
        if (n_lost < n):
            to_ret.append(tmp)
    return to_ret


def cache_images(cache_dir,func,**kw):
    """
    either caches or re-func to get an image transformaiton

    Args:
        cache_dir: where the cache is, or will be created
        func: functor (no arguments), re-reads everything if needed
        **kw: passed to CheckpointUtilities.multi_load
    Returns:
        cached imags
    """
    return CheckpointUtilities.multi_load(cache_dir,load_func=func,
                                          name_func=FEC_Util.name_func,**kw)

def get_coordinate_path(coords):    
    """
    Returns a path of 'coords' (assumed a single <x,y> coordinates of a 
    skeleton) such that:

    1) the start (endpoint) has the highest, second-lowest distance (the lowest
    distance is always +/- 1 pixel; the second lowest will be the greatest
    for an endpoint)
    
    2) all other points are separated from each other by at most 1 pixel

    Args:
        coords: the two-column array of pixel distances
    Returns:
        the sorted coordinate list
    """
    n_coords = len(coords)
    distances = scipy.spatial.distance_matrix(coords,coords,p=2)
    for i in range(n_coords):
        distances[i,i] = np.inf
    # check that the skeletonization is OK
    maximum_of_minimum_distances = np.sqrt(2)
    max_of_min = max(np.min(distances,axis=0))
    assert abs(max_of_min - maximum_of_minimum_distances) < 1e-6 , \
        "Skeletonization failed?"
    # POST: distances okay; all pixels at most 1 away in x and y
    # Now we need to decide on (possible arbitrary) endpoints. These should
    # be the two nodes with the largest *second* lowest distances (all have
    # at least one neighbor which is +/- 1 pixel; 'interior' nodes have at 
    # least two '1-pixel' neighbords
    second_lowest_distances = [sorted(row)[1] for row in distances]
    # sorted from low to high; what we want is the highest, second lowest
    sort_idx_second_highest = np.argsort(second_lowest_distances)
    endpoint = sort_idx_second_highest[-1]
    # POST: have endpoint. Add all the points with their two closest to the 
    # graph (except the endpoint, where we only add its closest)
    # create a graph of all the pixels
    G = nx.Graph()
    for i in range(n_coords):
        closest_nodes = np.argsort(distances[i])
        # add the closest
        G.add_edge(i,closest_nodes[0])
        if (i != endpoint):
            # also add the second closest for all 'interior' nodes...
            G.add_edge(i,closest_nodes[1])
    path = np.array(list(nx.dfs_preorder_nodes(G,endpoint)))
    return coords[path]

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    in_dir = "./data/"
    cache_dir_fmt = "./cache_{:s}/"
    cache_dir_raw = cache_dir_fmt.format("raw")
    out_dir = in_dir
    force_def = dict(force=True)
    GenUtilities.ensureDirExists(out_dir)
    images = ImageUtil.cache_images_in_directory(pxp_dir=in_dir,
                                                 cache_dir=cache_dir_raw,
                                                 limit=4)
    corrected_dir = cache_dir_fmt.format("corrected")
    label_dict = dict(cmap=plt.cm.spectral)
    transforms = [transform("corrected",correct_background),
                  transform("gaussian",blur_images),
                  transform("threshold",threshold_images),
                  transform("binarize",binarize_images),
                  transform("skeletonize",skeletonize_images),
                  transform("label",label_images,label_dict),
                  transform("skeleton_filter",skeleton_filter,label_dict)]
    last = images
    all_transforms = [last]
    for tx in transforms:
        tmp_dir = cache_dir_fmt.format(tx.name)
        last = cache_images(tmp_dir,func = lambda: tx.function(last),
                            **force_def)
        all_transforms.append(last)
    # fit a spline to the original data using each of the connected
    # regions from the skeletonization 
    image,skeleton = images[-1],last[-1]
    regions = measure.regionprops(skeleton.height)
    region = regions[0]
    coords = get_coordinate_path(region.coords)
    endpoint_coord = coords[0]
    n_coords = len(coords)
    x = coords[:,0]
    y = coords[:,1]
    idx = np.arange(n_coords)
    f_x = scipy.interpolate.interp1d(x=idx,y=x)
    f_y = scipy.interpolate.interp1d(x=idx,y=y)
    interp_idx = np.linspace(0,n_coords-1,num=10*n_coords,endpoint=True)
    x_interp = f_x(interp_idx)
    y_interp = f_y(interp_idx)
    plt.subplot(2,1,1)
    plt.imshow(image.height,origin='lower')
    plt.subplot(2,1,2)
    plt.plot(coords[:,1],coords[:,0],',')
    plt.plot(endpoint_coord[1],endpoint_coord[0],'go')
    plt.plot(coords[:,1],coords[:,0],'g-',alpha=0.3)
    plt.show()
    last_dir = tmp_dir
    # subtract the linear backround from each, save to a new cache 
    for i in range(len(images)):
        pipeline = [x[i] for x in all_transforms]
        # add in the keywords for the first image...
        imshow_kw = [dict(cmap=plt.cm.afmhot)] + \
                    [tx.imshow_kw for tx in transforms]
        fig = plot_with_background_corrected(pipeline,imshow_kw_list=imshow_kw)
        img_name = FEC_Util.name_func(i,images[i])
        out_name = "{:s}{:s}.png".format(last_dir,img_name)
        PlotUtilities.savefig(fig,out_name)

    
    
                                                 
if __name__ == "__main__":
    run()
