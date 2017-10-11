# force floating point division. Can still use integer with //
from __future__ import division
# other good compatibility recquirements for python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../../../../../../../../../")
from GeneralUtil.python import PlotUtilities,CheckpointUtilities
from GeneralUtil.python.Plot import Scalebar,Inset
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass,WeierstrassUtil
from GeneralUtil.python.IgorUtil import SavitskyFilter
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def run():
    """
    """
    landscape = CheckpointUtilities.lazy_load("./example_landscape.pkl")
    # make the landscape relative
    landscape.offset_energy(min(landscape.G_0))
    landscape.offset_extension(min(landscape.q))
    # get the landscape, A_z in kT. Note that we convert z->q, so it is
    # really A(q=z-A'/k)
    A_q_orig = landscape.A_z
    kT_to_kcal_mol = 0.593
    A_q_kT_orig = (A_q_orig * landscape.beta * kT_to_kcal_mol)
    n_points = 500
    A_q = SavitskyFilter(A_q_orig,n_points)
    A_q_kT = SavitskyFilter(A_q_kT_orig,n_points)
    q = SavitskyFilter(landscape.q,n_points)
    # filter the landscapes 
    # numerically differentiate
    to_y = lambda x: x * 1e12
    landscape_deriv_plot = to_y(np.gradient(A_q))
    # compare with the A' term. XXX should just save it...
    weighted_deriv_plot = to_y(landscape.A_z_dot)
    weighted_dA = to_y(landscape.A_z_dot * np.gradient(q))
    x_plot = q * 1e9
    energy_units = "kcal/mol"
    label_A_q_dot = r"$\dot{A}$"
    label_finite = "$dA$ from finite difference"
    label_work = r"$dA$ =<<$F$>> $dq$"
    label_dot = r"$\dot{A}$ =<<$F$>>"
    kw_weighted = dict(color='m')
    fig = PlotUtilities.figure((3.5,5))
    # # plot just A(q)
    ax_A_q = plt.subplot(3,1,1)
    color_energy = 'c'
    plt.plot(x_plot,A_q_kT,color=color_energy,label="$A$")
    PlotUtilities.lazyLabel("","$A$($q$) ("+ energy_units + ")","",
                            loc=(0.5,0.8),frameon=True)
    PlotUtilities.set_legend_kwargs(ax=ax_A_q,background_color='w',linewidth=0)
    PlotUtilities.no_x_label(ax_A_q)
    x0 = 14.5
    dx = 0.010
    xlim = [x0,x0+dx]
    zoom_x,zoom_y,ylim = Inset.slice_by_x(x_plot,A_q_kT,xlim)
    # add in some extra space for the scalebar 
    ylim_fudge = 0.1
    dy = ylim[1] - ylim[0]
    ylim = [ylim[0],ylim[1] + (ylim_fudge * dy)]
    lazy_common = dict(title_kwargs=dict(loc='left'))
    plt.plot(zoom_x,zoom_y,color=color_energy)
    # plot a zoomed in axis, to clarify why it probably goes wrong 
    axins = Inset.zoomed_axis(ax=ax_A_q,zoom=2000,xlim=xlim,ylim=ylim,
                              remove_ticks=True)
    mark_inset(parent_axes=ax_A_q,inset_axes=axins,loc1=2, loc2=3,
               ec="r")
    axins.plot(x_plot, A_q_kT,linewidth=1,color=color_energy)    
    # add in a scale bar for the inset. x goes from nm to pm (factor of 1000)
    unit_kw_x = dict(fmt="{:.0f}",value_function=lambda x: x*1000)
    common = dict(line_kwargs=dict(linewidth=1.0,color='k'))
    # round to 1 sig fig
    x_width = Scalebar.round_to_n_sig_figs(dx/2,1)
    y_width = Scalebar.round_to_n_sig_figs(dy * 0.68,1)
    unit_kw_x=dict(value_function=(lambda x: x*1000))
    Scalebar.scale_bar_rectangle_x(ax=axins,x_rel=0.5,y_rel=1.2,unit="pm",
                                   width=x_width,height_rel=0.3,
                                   unit_kwargs=unit_kw_x)
    Scalebar.scale_bar_rectangle_y(ax=axins,x_rel=-0.13,y_rel=0.53,
                                   unit="\ncal/mol",
                                   height=y_width,width_rel=0.2,font_color='w',
                                   unit_kwargs=unit_kw_x)
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    # # plot A_z_dot 
    ax_deriv_both = plt.subplot(3,1,2)
    # divide to get uN
    plt.plot(x_plot, landscape_deriv_plot * 1e9,color='k',
             label=label_finite)
    plt.plot(x_plot, weighted_dA * 1e9,label=label_work,**kw_weighted)
    PlotUtilities.lazyLabel("",
                            "$dA(q)$ ($\mathrm{\mu}$cal/mol)",
                            "$\Downarrow$ Calculate derivative",
                            **lazy_common)
    PlotUtilities.no_x_label(ax_deriv_both)
    # # plot A_z_dot, but just the weighted method (ie: not super wacky)
    ax_deriv_weighted = plt.subplot(3,1,3)
    plt.plot(x_plot,weighted_deriv_plot,linewidth=1,label=label_dot,
             **kw_weighted)
    title_last = "$\Downarrow$ Work-weighted method"
    PlotUtilities.lazyLabel("Extension (nm)","$\dot{A}(q)$ (pN)",
                            title_last,**lazy_common)
    PlotUtilities.label_tom(fig,axis_func=lambda ax: [ax[0]] + ax[2:],
                            loc=(-0.1,1.05),labels=PlotUtilities._lowercase)
    PlotUtilities.savefig(fig,"./finite_differences.png",
                          subplots_adjust=dict(hspace=0.2))

if __name__ == "__main__":
    run()
