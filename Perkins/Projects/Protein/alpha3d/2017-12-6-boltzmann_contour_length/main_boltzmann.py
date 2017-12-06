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

sys.path.append("../../../../../../")
from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from FitUtil.WormLikeChain.Python.Code import WLC_Utils
from scipy.interpolate import interp2d

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    _,data = FEC_Util.read_and_cache_pxp(directory="./",force=False)
    # offset and flip everything. All units are SI 
    x_offset = -4.62e-7
    f_offset = -1e-12
    max_ext = 60e-9
    for d in data:
        d.Force *= -1
        d.Force -= f_offset
        d.Separation -= x_offset
    # make a grid to convert to contour length space.
    desired_points_per_nm = 5
    ext = np.linspace(0,max_ext,num=max_ext * 1e9)
    min_L0 = 20e-9
    max_L0 = 55e-9
    num_points_for_resolution = (max_L0-min_L0) * 1e9 * desired_points_per_nm
    # create a grid for interpolating from extension and force to contour length
    L0 = np.linspace(min_L0,max_L0,num=num_points_for_resolution)
    Lp = 0.36e-9
    kT = 4.1e-21
    f_arr = []
    L0_arr = []
    ext_arr = []
    cutoff = 0.8
    print(L0,ext)
    for L0_tmp in L0:
        f = WLC_Utils.WlcNonExtensible(ext=ext,kbT=kT,Lp=Lp,L0=L0_tmp)
        valid_idx = np.where( ext < cutoff * L0_tmp)
        valid_f = f[valid_idx]
        f_arr.extend(valid_f)
        L0_arr.extend([L0_tmp for _ in valid_f])
        ext_arr.extend(ext[valid_idx])
    # make ext 2d...
    ext_flat = np.array(ext_arr)
    f_flat = np.array(f_arr)
    L0_flat = np.array(L0_arr)
    plt.plot(ext_flat,f_flat,'r.')
    plt.show()
    # POST: both f_arr and ext_2d are two dimensional. want to interpolate
    # from an ext/f pair to a contour length.
    key = data[0]
    print(L0_flat.shape,f_flat.shape,ext_flat.shape)
    interp = interp2d(x=ext_flat,y=f_flat,z=L0_flat,kind='linear',copy=True,
                      bounds_error=False)
    print(np.isfinite(key.Force).all())
    idx_i = np.where(key.Separation > 0)[0][0]
    idx_f = np.where(key.Separation > max_L0)[0]
    if (idx_f.size > 0):
        idx_f = idx_f[0]
    else:
        idx_f = -1
    print(idx_i,idx_f)
    m_slice = slice(idx_i,idx_f,1)
    L0_over_time = interp(key.Separation[m_slice],key.Force[m_slice])
    print(L0_over_time.shape)
    plt.subplot(2,1,1)
    plt.plot(key.Separation[m_slice],key.Force[m_slice])
    plt.subplot(2,1,2)
    plt.plot(key.Separation[m_slice],L0_over_time[0,:])
    plt.show()


if __name__ == "__main__":
    run()
