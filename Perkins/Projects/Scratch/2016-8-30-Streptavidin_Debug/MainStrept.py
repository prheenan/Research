# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../../../../")
from GeneralUtil.python import PlotUtilities



def run():
    """
    Makes plot of the streptavidin concentrations before and after
    aggresive spinning (14000rpm, 20 minutes)
    """
    # get the pre and post A280, A260/280, and concentration, if present 
    Pre = [ [1.482,0.57],
            [0.608,0.59],
            [1.430,0.56],
            [0.625,0.56],
            [1.431,0.57]]
    # post spinnng A280, A260/280
    Post = [ [1.034,0.58],
             [0.624,0.63],
             [1.479,0.60],
             [0.602,0.59],
             [1.398,0.57]]
    # see:
    # 1.Biotin binding proteins.
    #Available at:
    #abcam.com/index.html?pageconfig=resource&rid=12723.
    #(Accessed: 30th August 2016)
    ExtinctionPerMoleCm = 41326
    MolecularMassDaltons = 52800
    # nanodrop give 10mm effective absorbances
    PathLengthCm = 1
    # coefficienct to converty from intensity to mg/ml, using pp 133 of
    # notebook #2 (Beers law)
    AbsorbanceToCocnentrationMolesPerLiter = \
        1/(ExtinctionPerMoleCm*PathLengthCm)
    LitersToMl = 1000
    DaltonToKg = 1.66e-27
    MolesToMolecules = 6.02e23
    KbToMg = 1e6
    MolesToMg = MolesToMolecules * MolecularMassDaltons * DaltonToKg * KbToMg
    print(MolesToMg)
    MolesPerLiterToMgPerMl = MolesToMg/LitersToMl
    coeff = AbsorbanceToCocnentrationMolesPerLiter*MolesPerLiterToMgPerMl
    GetConc = lambda arr,skip: np.array([a[0] * coeff for a in arr[skip:]])
    # skip the first data point, which was just a control
    skip = 1
    PreConc = GetConc(Pre,skip)
    PostConc = GetConc(Post,skip)
    fig = PlotUtilities.figure()
    MinV = min(min(PreConc),min(PostConc))
    MaxV = max(max(PreConc),max(PostConc))
    Lim = [MinV*0.95,MaxV*1.05]
    # error from
    # Thermo Scientific. NanoDrop 1000 Spectrophotometer. Page 38
    # Available at: nanodrop.com/library/nd-1000-v3.7-users-manual-8.5x11.pdf.
    # (Accessed: 30th August 2016)
    plt.errorbar(PreConc,PostConc,yerr=0.1,fmt='ro',label="Data")
    plt.plot(Lim,Lim,
             color='b',linewidth=2,alpha=0.5,label="before=after")
    plt.xlim(Lim)
    plt.ylim(Lim)
    PlotUtilities.lazyLabel("Concentration before spinning [mg/mL]",
                            "Concentration after spinning [mg/mL]",\
                            "[Protein] invariant, post-20', 14000 RPM spin",
                            frameon=True)
    PlotUtilities.savefig(fig,"Conc.png")


if __name__ == "__main__":
    run()
