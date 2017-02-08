# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys


def add_stretch(x,f,start,end,spring):
    f_copy = f.copy()
    idx_stretch = np.where( (x > start) & 
                            (x < end))
    x_stretch = x[idx_stretch]
    f_stretch = (x_stretch - x_stretch[0])**2 * spring
    return idx_stretch,f_stretch

def make_force_extension_curve(x,array_of_stretch_kwargs,DecayConst,snr):
    f = (1-np.exp(-x/DecayConst))
    for kwargs in array_of_stretch_kwargs:
        idx,f_at_idx = add_stretch(x=x,f=f,**kwargs)
        f[idx] += f_at_idx
    # add in uniform noise
    noise_ampltude = np.sqrt(1/snr)
    f += (np.random.normal(size=f.size)-0.5)* 2 * noise_ampltude
    return f

def normalized_force_extension(max_x=1,decay_const=1/50,points=1024,
                               spring_const=10,snr=100):
    x,dx = np.linspace(0,max_x,points,retstep=True)
    SpringStretch = spring_const
    DecayConst = decay_const
    stretch_kwargs = [
        # adhesion peak...
        dict(start=0.04,end=0.1,spring=SpringStretch*50),
        # next, add in something that looks (vaguely) like a WLC strech
        dict(start=0.2,end=0.4,spring=SpringStretch),
        dict(start=0.35,end=0.6,spring=SpringStretch*1.2),
        dict(start=0.55,end=0.9,spring=SpringStretch/2)]
    # get the force
    f = make_force_extension_curve(x,stretch_kwargs,DecayConst,snr)
    # convert to 0-1 (ish)
    f -= np.median(f)
    f /= np.max(f)
    # get the maximum force in 
    return x,dx,stretch_kwargs,f
