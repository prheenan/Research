# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,copy

sys.path.append("../../../../../../../../../../../")
sys.path.append("../../../")
from Util import IoUtil,AnalysisClasses
from Research.Perkins.AnalysisUtil.Images import PolymerTracing,ImageUtil
from GeneralUtil.python import PlotUtilities,CheckpointUtilities,GenUtilities
from GeneralUtil.python.Plot import Scalebar
from skimage.morphology import square,dilation
from matplotlib.colors import Normalize

def _trace_x_y_plot(t,um_per_px):
    x,y = t.inf.x_y_abs
    x_plot,y_plot = x * um_per_px, y * um_per_px
    return x_plot,y_plot

def plot_trace(ax,t,um_per_px):
    if (t.inf is not None):
        x_plot,y_plot = _trace_x_y_plot(t,um_per_px)         
        ax.plot(x_plot,y_plot,linewidth=0.25)

def um_per_px(i):
    m_per_px = i.m_per_px
    um_per_px = m_per_px * 1e6
    return um_per_px

def overlay_traces(fig,i,output_path):
    im = i.image
    ax = plt.subplot(1,1,1)
    ImageUtil.image_plot(im,ax=ax,fig=fig)
    for t in i.worm_objects:
        plot_trace(ax,t,um_per_px(i))

def scalebar(ax,square_um,color):
    scalebar_width_um = square_um 
    unit_kwargs = dict(value_function=lambda x: 1000*x)
    x_font,_ = Scalebar.font_kwargs_modified(x_kwargs=dict(color=color))
    line_kwargs = dict(**Scalebar.def_line_kwargs)
    line_kwargs['color'] = color
    Scalebar.x_scale_bar_and_ticks_relative(unit="nm",width=scalebar_width_um,
                                            offset_x=0.5,offset_y=0.95,
                                            unit_kwargs=unit_kwargs,
                                            font_kwargs=x_font,
                                            ax=ax,
                                            line_kwargs=line_kwargs)
    PlotUtilities.no_x_anything(ax)
    PlotUtilities.no_y_anything(ax)
    ax.axis('off')



def single_trace(fig,i,trace,add_trace=True,add_scalebar=True,
                 skeletonize=4,**kw):
    ax = plt.subplot(1,1,1)
    um_per_px_i = um_per_px(i)
    image = copy.deepcopy(i.image)
    # copy the height, in case we modify it...
    image.height = copy.deepcopy(image.height)
    skeleton_bool = skeletonize is not None and skeletonize > 0
    if (skeleton_bool):
        # get the pixel representation...
        x_px,y_px = _trace_x_y_plot(trace,1)
        skeleton = np.zeros(image.height.shape)
        for x_tmp,y_tmp in zip(x_px,y_px):
            skeleton[int(x_tmp),int(y_tmp)] = 1
        skeleton = dilation(skeleton, square(skeletonize))
        cmap = kw['cmap']
        data = copy.deepcopy(image.height)
        data *= skeleton
        min_d,max_d = np.min(data),np.max(data)
        data -= min_d
        data /= (max_d-min_d)
        x_thresh,y_thresh = np.where(skeleton < 1)
        rgba_img = cmap(data)
        rgba_img[x_thresh,y_thresh,:-1] = 1
        rgba_img[x_thresh,y_thresh,-1] = 0
        plt.imshow(rgba_img,extent=[0,2,0,2],interpolation='bicubic',
                   **kw)
    else:
        image_plot = ImageUtil.image_plot(image,ax=ax,fig=fig,imshow_kwargs=kw)
    square_um = 0.30
    x,y = _trace_x_y_plot(trace,um_per_px_i)
    x0 = np.mean(x)
    y0 = np.mean(y)
    xlim = [x0-square_um,x0+square_um]
    ylim = [y0-square_um,y0+square_um]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if (add_scalebar):
        scalebar(ax,square_um,color=('k' if not skeleton_bool else 'k'))
    if (add_trace):
        plot_trace(ax,trace,um_per_px_i)

    
def plot_single_trace(fig,image,trace,add_scalebar=False,**kw):
    single_trace(fig,image,trace,
                 add_trace=False,
                 add_scalebar=add_scalebar,**kw)

def plot_image(image,output_dir):
    name = GenUtilities.file_name_from_path(image.image_path)
    output_path = output_dir + name + "all_traces.png"
    # make a plot with everything
    fig = PlotUtilities.figure()
    overlay_traces(fig,image,output_path)
    PlotUtilities.savefig(fig,output_path)
    # make a plot with each image
    imshow_kwargs = dict(**ImageUtil.def_imshow_kw)
    imshow_kwargs['cmap'] = plt.cm.afmhot_r
    for j,trace in enumerate(image.worm_objects):
        if (trace.inf is None):
            continue
        fig = PlotUtilities.figure(frameon=False)
        plot_single_trace(fig,image,trace,skeletonize=False,
                          add_scalebar=True,**imshow_kwargs)
        output_name = output_path[:-3] + "{:d}_.pdf".format(j)
        PlotUtilities.savefig(fig,output_name,transparent=True)

def run(in_dir):
    """
    Args:
        in_dir: the input directory to operate on.  
    """
    cache_dir = IoUtil._traces_dir(in_dir)
    output_dir = "./" #IoUtil._plot_dir(in_dir)
    images = CheckpointUtilities.lazy_multi_load(cache_dir)
    for i,im in enumerate(images):
        plot_image(im,output_dir + "{:d}_".format(i))
        exit
    pairs = [[12,2,"part of e"],
             [22,2,"l"],
             [2,0,"part of a"],
             [15,1,"part of a"],
             [2,1,"part of d,j"],
             [2,5,"part of y"],
             [11,2,"part of n"],
             [15,1,"part of n"],
             [2,5,"v"],
             [15,3,"i"],
             [17,1,"i"],
             [4,0,"s"],
             [16,1,"s"],
             [3,1,"o"],
             [16,2,"o"],
             [8,3,"part of r"],
             [11,2,"part of r"],
             [13,5,"part of r"],
             [18,0,"butterfly"]]
    images_traces = [ [images[i],images[i].worm_objects[j]] for i,j,_ in pairs]
    imshow_kwargs = dict(**ImageUtil.def_imshow_kw)
    imshow_kwargs['cmap'] = plt.cm.afmhot
    for i,(im,t) in enumerate(images_traces):
        fig = PlotUtilities.figure(figsize=(5,5),frameon=False)
        plot_single_trace(fig,im,t,skeletonize=None,**imshow_kwargs)
        output_name = "./mom/" + "{:d}_{:s}.pdf".format(i,pairs[i][-1])
        PlotUtilities.savefig(fig,output_name,transparent=True)   
    # make a single one with a scale bar, for reference. 
    fig = PlotUtilities.figure(figsize=(5,5),frameon=False)
    plot_single_trace(fig,*images_traces[0],
                      add_scalebar=True,**imshow_kwargs)
    
        


if __name__ == "__main__":
    run(IoUtil.get_directory_command_line())
