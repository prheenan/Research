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
from scipy.interpolate import splprep, splev, interp1d,UnivariateSpline

from sklearn.neighbors import NearestNeighbors
from route.postman import single_chinese_postman_path

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
    n_neightbors = 2
    # sort the data so the endpoint is first?
    print(endpoint)
    sorted_idx = list(np.arange(endpoint,n_coords)) + \
                 list(np.arange(0,endpoint))
    sorted_idx= np.array(sorted_idx)
    distances = distances[sorted_idx]
    coords = coords[sorted_idx]
    for i in range(n_coords):
        dist_tmp = distances[i]
        closest_nodes = np.argsort(dist_tmp)
        # add the closest N
        dist_sort = dist_tmp[closest_nodes]
        G.add_edge(i,closest_nodes[0],weight=dist_sort[0])
        if (i != endpoint):
            G.add_edge(i,closest_nodes[1],weight=dist_sort[1])
    graph,path = single_chinese_postman_path(G)
    print(path,n_coords)
    for i in range(len(path)):
        print(len(set(path[:i])),i,n_coords)
    path = path[:152]
    """
    see: 
https://stackoverflow.com/questions/18794308/algorithm-to-cover-all-edges-given-starting-node

https://networkx.github.io/documentation/networkx-1.9.1/reference/generated/networkx.algorithms.matching.max_weight_matching.html#networkx.algorithms.matching.max_weight_matching

    also https://groups.google.com/forum/#!topic/networkx-discuss/NxbsY2dzkNk
    
    https://healthyalgorithms.com/2009/03/23/aco-in-python-minimum-weight-perfect-matchings-aka-matching-algorithms-and-reproductive-health-part-4/
    """
    coords_x = np.array(coords[:,0])
    coords_y = np.array(coords[:,1])

    return coords[path]

def snake_fit(image,initial,w_line=5,w_edge=0,max_px_move=1,beta=1,gamma=0.1):
    to_fit = image
    min_x,max_x = np.min(initial[:,0]),np.max(initial[:,0])
    min_y,max_y = np.min(initial[:,1]),np.max(initial[:,1])
    fudge_x = int(np.ceil((max_x-min_x) * 0.1))
    fudge_y = int(np.ceil((max_y-min_y) * 0.1))
    lower_x = 0#max(0,min_x-fudge_x)
    lower_y = 0#max(0,min_y-fudge_y)
    #to_fit = to_fit[lower_x:max_x+fudge_x,
    #                lower_y:max_y+fudge_y]
    initial_x_shifted = initial[:,0]-lower_x
    initial_y_shifted = initial[:,1]-lower_y
    initial = np.array((initial_x_shifted,initial_y_shifted)).T
    min_image,max_image = np.min(to_fit),np.max(to_fit)
    to_fit =  ((to_fit - min_image)/(max_image - min_image)) * 256
    to_fit = to_fit.astype(np.uint8)
    initial_snake = initial.astype(np.float64)
    snake = active_contour(to_fit,convergence=1e-3,max_iterations=5e3,
                           snake=initial_snake,w_line=w_line,
                           w_edge=w_edge,beta=beta,gamma=gamma,
                           bc='fixed',max_px_move=max_px_move)
    snake[:,0] += lower_x
    snake[:,1] += lower_y
    return snake


def plot_fitting(image,coords,snake_coords=None):
    endpoint_coord = coords[0]
    n_coords = len(coords)
    x = coords[:,0]
    y = coords[:,1]
    idx = np.arange(n_coords)
    plt.imshow(image.T,origin='lower')
    plt.plot(coords[:,0],coords[:,1],',')
    plt.plot(endpoint_coord[0],endpoint_coord[1],'go')
    plt.plot(coords[:,0],coords[:,1],'r-',alpha=0.3)
    plt.xlim(min(coords[:,0])*0.8,max(coords[:,0]*1.1))
    plt.ylim(min(coords[:,1])*0.8,max(coords[:,1]*1.1))
    if (snake_coords is not None):
        plt.plot(snake_coords[:,0],snake_coords[:,1],'r.-',linewidth=0.3)

def get_spline_obj(image,coords,fudge=3,k=3,smooth_f_n=2):
    """
    Returns:
        tuple of 
    """
    n_coords = len(coords)
    coords_x = coords[:,0]
    coords_y = coords[:,1]
    # threshold the image outside of the skeleton area of interest
    image_thresh = copy.deepcopy(image)
    zero_x_low = coords_x - fudge
    zero_x_high = coords_x + fudge
    zero_y_low = coords_y - fudge
    zero_y_high = coords_y + fudge
    # get a mask with ones in the region...
    m_arr = np.zeros(image.height.shape)
    for x_l,x_h,y_l,y_h in zip(zero_x_low,zero_x_high,
                               zero_y_low,zero_y_high):
        m_arr[x_l:x_h,y_l:y_h] = 1
    image_thresh.height *= m_arr
    where_non_zero_image_xy =np.where(image_thresh.height > 0)
    where_non_zero_x,where_non_zero_y = where_non_zero_image_xy
    """
    # Essentially trying to parameterize a curve basd on the skeleton. see: 
stackoverflow.com/questions/31464345/fitting-a-closed-curve-to-a-set-of-points
    also:
stackoverflow.com/questions/32046582/spline-with-constraints-at-border/32421626#32421626
stackoverflow.com/questions/36830942/reordering-image-skeleton-coordinates-to-make-interp1d-work-better
    https://stackoverflow.com/questions/41659075/how-to-specify-the-number-of-knot-points-when-using-scipy-splprep

    ... also ...

    https://stackoverflow.com/search?q=parametric+image+fit

    especially:

    stackoverflow.com/questions/22556381/approximating-data-with-a-multi-segment-cubic-bezier-curve-and-a-distance-as-wel/22582447#22582447

    and

    stackoverflow.com/questions/22556381/approximating-data-with-a-multi-segment-cubic-bezier-curve-and-a-distance-as-wel/22582447#22582447
    """
    n_non_zero = where_non_zero_x.size
    # get the projection of the data onto the skeleton
    skel_idx = []
    weights = []
    for i,(x,y) in enumerate(zip(where_non_zero_x,where_non_zero_y)):
        diff = np.sqrt((coords_x-x)**2 + (coords_y-y)**2)
        closest_skeleton_idx = np.argmin(diff)
        skel_idx.append(closest_skeleton_idx)
        weights.append(image_thresh.height[x,y])
    # sort the projection array by the distance along...
    sort_idx = np.argsort(skel_idx)
    sorted_image_x = where_non_zero_x[sort_idx]
    sorted_image_y = where_non_zero_y[sort_idx]
    f_excess = where_non_zero_y.size/len(coords_x)
    # weight; sum by default is such that sum is M (being the size). use that
    weights = np.array(weights)[sort_idx]
    weights /= np.sum(weights)
    weights *= weights.size
    n_points = sorted_image_x.size
    tck,_ = splprep([sorted_image_x,sorted_image_y],w=weights,per=0,
                    s=f_excess*n_points*smooth_f_n,k=k,u=None,quiet=0)
    return image_thresh,tck

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
    last_dir = tmp_dir
    # fit a spline to the original data using each of the connected
    # regions from the skeletonization 
    for image,skeleton in zip([images[-2]],[last[-2]]):
        regions = measure.regionprops(skeleton.height)
        region = regions[0]
        coords = get_coordinate_path(region.coords)
        unew = np.linspace(0,1,endpoint=True,num=10*max(coords.shape))
        image_thresh,tck = get_spline_obj(image,coords)
        out = splev(unew,tck)
        out_x ,out_y = out[0],out[1]
        min_x,max_x = min(out_x),max(out_x)
        range_v = (max_x-min_x) * 0.2
        xlim = lambda : plt.xlim([min_x-range_v,max_x+range_v])
        ylim = lambda : plt.ylim([min(out_y)-range_v,max(out_y+range_v)])
        plt.subplot(2,1,1)
        plt.imshow(skeleton.height.T)
        plt.plot(coords[:,0],coords[:,1],'r-')
        xlim()
        ylim()
        plt.subplot(2,1,2)
        plt.imshow(image_thresh.height.T)
        plt.plot(out_x,out_y,'r-')
        plt.plot(coords[0,0],coords[0,1],'go')
        xlim()
        ylim()
        plt.show()
    fig = PlotUtilities.figure()
    plot_fitting(snake_input.height,coords)
    plt.plot(interp_x(ind_sort), interp_y(ind_sort), 'r,')
    out_name = "{:s}_fit.png".format(out_dir)
    PlotUtilities.savefig(fig,out_name)
    exit(1)
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
