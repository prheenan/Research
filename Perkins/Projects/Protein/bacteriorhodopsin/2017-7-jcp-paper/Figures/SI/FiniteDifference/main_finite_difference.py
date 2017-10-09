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
    n_points = 2000
    A_q = SavitskyFilter(A_q_orig,n_points)
    A_q_kT = SavitskyFilter(A_q_kT_orig,n_points)
    q = SavitskyFilter(landscape.q,n_points)
    # filter the landscapes 
    # numerically differentiate
    to_y = lambda x: x * 1e12
    landscape_deriv_plot = to_y(np.gradient(A_q)/np.gradient(q))
    # compare with the A' term. XXX should just save it...
    weighted_deriv_plot = to_y(landscape.A_z_dot)
    x_plot = q * 1e9
    energy_units = "kcal/mol"
    label_A_q_dot = r"$\dot{A}$"
    label_finite = label_A_q_dot + r" from finite difference"
    label_work = r"{:s}$ =<<F>>$".format(label_A_q_dot)
    kw_weighted = dict(color='m',label=label_work)
    fig = PlotUtilities.figure((3.5,5))
    # # plot just A(q)
    ax_A_q = plt.subplot(3,1,1)
    color_energy = 'c'
    plt.plot(x_plot,A_q_kT,color=color_energy,label="$A$")
    PlotUtilities.lazyLabel("","Helmholtz A ("+ energy_units + ")","",
                            loc=(0.5,0.8),frameon=True)
    PlotUtilities.set_legend_kwargs(ax=ax_A_q,background_color='w',linewidth=0)
    PlotUtilities.no_x_label(ax_A_q)
    x0 = 14.5
    dx = 0.05
    xlim = [x0,x0+dx]
    zoom_x,zoom_y,ylim = Inset.slice_by_x(x_plot,A_q_kT,xlim)
    # add in some extra space for the scalebar 
    ylim_fudge = 0.7
    dy = ylim[1] - ylim[0]
    ylim = [ylim[0],ylim[1] + (ylim_fudge * dy)]
    lazy_common = dict(title_kwargs=dict(loc='left'))
    plt.axvspan(*xlim,color=color_energy,alpha=0.3,edgecolor="None")
    plt.plot(zoom_x,zoom_y,color='r')
    # plot a zoomed in axis, to clarify why it probably goes wrong 
    axins = Inset.zoomed_axis(ax=ax_A_q,zoom=500,xlim=xlim,ylim=ylim,
                              remove_ticks=False)
    axins.plot(x_plot, A_q_kT,linewidth=0.1,color='r')    
    # add in a scale bar for the inset. x goes from nm to pm (factor of 1000)
    unit_kw_x = dict(fmt="{:.0f}",value_function=lambda x: x*1000)
    common = dict(line_kwargs=dict(linewidth=1.0,color='k'))
    # round to ~10s of pm
    x_width = np.around(dx/3,2)
    y_width = np.around(dy/3,1)
    x_kw = dict(width=x_width,unit="pm",unit_kwargs=unit_kw_x,
                fudge_text_pct=dict(x=0.2,y=-0.2),**common)
    y_kw = dict(height=y_width,unit=energy_units,
                unit_kwargs=dict(fmt="{:.1f}"),**common)
    ylim = axins.get_ylim()
    y_diff = ylim[1] - ylim[0]
    Scalebar._scale_bar_rectangle(ax=axins,x=0.5,
                                  s="{:.0f} pm".format(x_width*1000),
                                  y=1.3,width=x_width,
                                  height=y_diff * 0.4)
    # # plot A_z_dot 
    ax_deriv_both = plt.subplot(3,1,2)
    # divide by 1000 to get uN
    plt.plot(x_plot,landscape_deriv_plot/1e3,color='k',
             label=label_finite)
    plt.plot(x_plot,weighted_deriv_plot/1e3,**kw_weighted)
    PlotUtilities.lazyLabel("",
                            "$\dot{A}(q)$ (nN)",
                            "$\Downarrow$ Determine derivative (both methods)",
                            **lazy_common)
    PlotUtilities.no_x_label(ax_deriv_both)
    # # plot A_z_dot, but just the weighted method (ie: not super wacky)
    ax_deriv_weighted = plt.subplot(3,1,3)
    plt.plot(x_plot,weighted_deriv_plot,linewidth=1,**kw_weighted)
    title_last = "$\Downarrow$ Work-weighted method is reasonable "
    PlotUtilities.lazyLabel("Extension (nm)","$\dot{A}(q)$ (pN)",
                            title_last,**lazy_common)
    PlotUtilities.label_tom(fig,axis_func=lambda ax: [ax[0]] + ax[2:] )        
    PlotUtilities.savefig(fig,"./finite_differences.png",
                          subplots_adjust=dict(hspace=0.2))

if __name__ == "__main__":
    run()
