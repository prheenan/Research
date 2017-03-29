# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os,string

sys.path.append("../../../../../../../../")
from GeneralUtil.python import PlotUtilities
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util,FEC_Plot
from Research.Personal.EventDetection.Util.InputOutput import \
    read_and_cache_file
from Research.Personal.EventDetection.Util import Analysis,Plotting
# /!\ note the 'SVG' function also in svgutils.compose
import matplotlib.gridspec as gridspec
from Research.Personal.EventDetection._2SplineEventDetector import Detector

def f_plot_y(y):
    return y*1e12

def plot(interp,split_fec,f,xlim_rel_start,xlim_rel_delta):
    time_sep_force = f(split_fec)
    x_plot,y_plot = Plotting.plot_format(time_sep_force)
    n_filter_points = 1000
    x_raw = time_sep_force.Time
    y_raw = time_sep_force.Force
    interp_raw = interp(x_raw)
    diff_raw = y_raw - interp_raw
    stdev = Analysis.local_stdev(diff_raw,n=split_fec.tau_num_points)
    xlims_rel = [ [i,i+xlim_rel_delta] for i in  xlim_rel_start]
    # convert to plotting units
    n = x_plot.size
    slices_abs = [ slice(int(i*n),int(f*n),1) for i,f in xlims_rel ]
    x_plot_slices = [ x_plot[s] for s in slices_abs ]
    diff_raw_slices = [diff_raw[s] for s in slices_abs]
    max_raw_diff_slices = [max(d) for d in diff_raw_slices]
    min_raw_diff_slices = [min(d) for d in diff_raw_slices]
    range_raw_diff_slices = np.array([min(min_raw_diff_slices),
                                      max(max_raw_diff_slices)])
    range_raw_diff = np.array([min(diff_raw),max(diff_raw)])
    range_plot_diff = f_plot_y(range_raw_diff*1.1)
    xlim_abs = [ [min(x),max(x)] for x in x_plot_slices]
    n_plots = len(x_plot_slices)
    # set up the plot styling
    style_approach = dict(color='b')
    style_raw = dict(alpha=0.3,**style_approach)
    style_interp = dict(linewidth=3,**style_approach)
    colors = ['r','m','k']
    style_regions = [] 
    for c in colors:
        style_tmp = dict(**style_raw)
        style_tmp['color'] = c
        style_regions.append(style_tmp)
    gs = gridspec.GridSpec(3,2*n_plots)
    plt.subplot(gs[0,:])
    plt.plot(x_plot,y_plot,**style_raw)
    plt.plot(x_plot,f_plot_y(interp_raw),**style_interp)
    PlotUtilities.lazyLabel("Time (s)","Force (pN)","")
    PlotUtilities.x_label_on_top()
    ax_diff = plt.subplot(gs[1,:])
    plt.plot(x_plot,f_plot_y(diff_raw),**style_raw)
    PlotUtilities.lazyLabel("","Residual (pN)","")
    PlotUtilities.no_x_label()
    # highlight all the residual regions in their colors
    for style_tmp,slice_tmp in zip(style_regions,slices_abs):
        plt.plot(x_plot[slice_tmp],f_plot_y(diff_raw)[slice_tmp],**style_tmp)
    plt.ylim(range_plot_diff)
    # plot all the subregions
    for i in range(n_plots):
        xlim_tmp = xlim_abs[i]
        ylim_tmp = range_plot_diff
        """
        plot the raw data
        """
        offset_idx = 2*i
        ax_tmp = plt.subplot(gs[-1,offset_idx])
        diff_tmp = diff_raw_slices[i]
        # convert everything to plotting units
        diff_plot_tmp =f_plot_y(diff_tmp) 
        x_plot_tmp = x_plot_slices[i]
        style_tmp = style_regions[i]
        plt.plot(x_plot_tmp,diff_plot_tmp,**style_tmp)
        PlotUtilities.no_x_anything()
        if (i != 0):
            PlotUtilities.no_y_label()
            PlotUtilities.tickAxisFont()
        else:
            PlotUtilities.lazyLabel("","Force (pN)","")
        plt.xlim(xlim_tmp)
        plt.ylim(ylim_tmp)
        PlotUtilities.zoom_effect01(ax_diff, ax_tmp, *xlim_tmp)
        if (i == 0):
            # make a scale bar for this plot
            time = 25e-3
            string = "{:d} ms".format(int(time*1000))
            PlotUtilities.scale_bar_x(np.mean(xlim_tmp),0.8*max(ylim_tmp),
                                      s=string,width=time,fontsize=15)
        """
        plot the histogram
        """
        plt.subplot(gs[-1,offset_idx+1])
        if (i == 0):
            PlotUtilities.xlabel("Count")
        color = style_tmp['color']
        # overlay a box plot on top
        n,_,_ = plt.hist(diff_plot_tmp,orientation="horizontal",**style_tmp)
        plt.boxplot(diff_plot_tmp,vert=True,positions=[np.max(1.2*n)],
                    manage_xticks=False,widths=(40),
                    medianprops=dict(linewidth=3,color='k'),
                    whiskerprops=dict(linewidth=1,**style_tmp),
                    flierprops=dict(marker='o',markersize=2,**style_tmp),
                    boxprops=dict(linewidth=3,fillstyle='full',**style_tmp))
        plt.ylim(ylim_tmp)
        plt.xlim([0,max(n)*1.6])
        PlotUtilities.tickAxisFont()
        PlotUtilities.no_y_label()

def run(base="./"):
    """
    
    """
    name = "examples.pdf"
    data_base = base + "data/"
    file_name = "fast_unfolding"
    kw = dict(cache_directory=data_base,force=False)
    file_path = data_base + file_name +".csv"
    fec = read_and_cache_file(file_path,**kw)
    split_fec = Analysis.zero_and_split_force_extension_curve(fec)
    xlim_rel_start_appr = [0.1,0.4,0.8]
    xlim_rel_start_retr = [0.2,0.4,0.54]
    xlim_rel_delta = 0.02
    interp = split_fec.approach_spline_interpolator()
    interp_retract = split_fec.retract_spline_interpolator()
    f_appr = lambda x: x.approach
    f_retr = lambda x: x.retract
    fig = PlotUtilities.figure((10,8))
    plot(interp_retract,split_fec,f=f_retr,xlim_rel_start=xlim_rel_start_retr,
         xlim_rel_delta=xlim_rel_delta)
    PlotUtilities.savefig(fig,"./out_retr.png")
    fig = PlotUtilities.figure((10,8))
    plot(interp,split_fec,f=f_appr,xlim_rel_start=xlim_rel_start_appr,
         xlim_rel_delta=xlim_rel_delta)
    PlotUtilities.savefig(fig,"./out.png")



if __name__ == "__main__":
    run()
