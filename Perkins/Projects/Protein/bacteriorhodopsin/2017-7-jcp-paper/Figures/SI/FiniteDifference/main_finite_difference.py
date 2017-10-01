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
from FitUtil.EnergyLandscapes.InverseWeierstrass.Python.Code import \
    InverseWeierstrass

def run():
    """
    """
    landscape = CheckpointUtilities.lazy_load("./example_landscape.pkl")
    # make the landscape relative
    landscape.offset_energy(min(landscape.G_0))
    landscape.offset_extension(min(landscape.q))
    # get the landscape, A_z in kT. Note that we convert z->q, so it is
    # really A(q=z-A'/k)
    A_q = landscape.A_z
    A_q_kT = (A_q * landscape.beta)
    # numerically differentiate
    to_y = lambda x: x * 1e12
    landscape_deriv_plot = to_y(np.gradient(A_q)/np.gradient(landscape.q))
    # compare with the A' term. XXX should just save it...
    weighted_deriv_plot = to_y(landscape.A_z_dot)
    x_plot = landscape.q * 1e9
    plt.subplot(3,1,1)
    plt.plot(x_plot,A_q_kT)
    plt.plot(x_plot,landscape.G_0 * landscape.beta)
    PlotUtilities.lazyLabel("Extension (nm)","A(q) (kT/nm)","")
    plt.subplot(3,1,2)
    # divide by 1000 to get uN
    plt.plot(x_plot,landscape_deriv_plot/1e6)
    plt.plot(x_plot,weighted_deriv_plot/1e6)
    PlotUtilities.lazyLabel("Extension (nm)",
                            "$\dot{A}(q)$ ($\mathrm{\mu}$N)","")
    plt.subplot(3,1,3)
    plt.plot(x_plot,weighted_deriv_plot)
    PlotUtilities.lazyLabel("Extension (nm)","$\dot{A}(q)$ (pN)","")
    plt.show()

if __name__ == "__main__":
    run()
