# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import Bio.SeqUtils
sys.path.append("../../../../../../../")
import GeneralUtil.python.PlotUtilities as pPlotUtil

def LinearRadiusOfGyration(a,L):
    """
    Equation 7 from:
    Latulippe, David R., and Andrew L. Zydney. 
    "Radius of Gyration of Plasmid DNA Isoforms from Static Light Scattering." 
    Biotechnology and Bioengineering 107, 2010
    """
    return a * np.sqrt(( L/(3*a) -1 + 2*a/L 
                         -2 * (a/L)**2 * (1 - np.exp(-L/a))))

def CircularRadiusOfGyration(a,L):
    """
    Equation 8 from:
    Latulippe, David R., and Andrew L. Zydney. 
    "Radius of Gyration of Plasmid DNA Isoforms from Static Light Scattering." 
    Biotechnology and Bioengineering 107, 2010
    """
    alpha = 0.5
    k2 = 0.5
    k3 = 0.25
    a_alph = a * alpha
    return np.sqrt( (2*a + (22*a**2)/(3*L)) * \
                    (L/12 - 2*(a_alph)**2/L +  8 * (a_alph)**3/(3*(L**2)))\
                    - a**2 + 4*alpha*(a**3)/L + 8*(a_alph**3)/L * \
                    (1/3 - alpha/6 + k2 * alpha**2/5 + k3*(alpha**3)/6))


def run():
    # get our actual sequence weight
    Primer = "AGAGTGGTCCTA"
    # get the sequence from 1607 to 3520, tack on the primer
    with open("mp13_plasmid_plasmid_seq.txt") as f:
        Seq = Primer + "".join([l for l in f.readlines()])[1607:3520]
    MolecularWeight = Bio.SeqUtils.molecular_weight(Seq,
                                                    seq_type="DNA",
                                                    double_stranded=True,
                                                    circular=True)
    LengthDoubleStrandedBasePairs = len(Seq)
    MassGramsPerMolecule= MolecularWeight/6.022e23
    MoleculesPerCircularMicrogram = 1e-6/MassGramsPerMolecule
    KbpPerConstr = len(Seq)/1000
    # #get the mean distance between molecules on an afm slide
    radius = 0.010
    Area = np.pi * (radius**2)
    # set up an array of DNA amounts loaded
    MicrogramsLoaded = np.logspace(-3,2,num=2000)
    # convert the micrograms into a concentration, assuming LoadVolume is the
    # volume of the DNA
    LoadVolumeMicroliters = 20
    LoadConcentrationNgPerUl = (MicrogramsLoaded * 1e3)/LoadVolumeMicroliters
    # get many molecules that mass translates into
    MoleculesLoaded = MoleculesPerCircularMicrogram * MicrogramsLoaded
    Avogadro = 6.022e23
    MolarityLoaded =(MoleculesLoaded/(LoadVolumeMicroliters * 1e-6))/Avogadro
    NanoMolarLoaded = MolarityLoaded * 1e9
    # get the expected mean distance between molecules, assuming 100%
    # binding efficiency
    MeanDist = np.sqrt(Area/MoleculesLoaded)
    # convert the distance to nanometers
    MeanDistNano = MeanDist * 1e9
    # get the DNA Size in nm
    NanometersPerBp = 0.338
    DNASizeNanoMeters = NanometersPerBp * KbpPerConstr * 1000
    # maxmium concentration fro QiaQuick is about 300ng/uL in 35uL for 10ug
    # anything higher than this is *not* recommended
    MaxConcNgPerUl = 300
    # the tip radius is O(10nm), persistence length is 43nm, consider that
    # as a 'length' of the polymer to find the radius of gyration
    Lp_nm = 43
    NumPersistence = DNASizeNanoMeters/Lp_nm
    # at least N radii of gyration away. Should be more than 2 to prevent
    # overlap
    RadiusOfGyr = LinearRadiusOfGyration(Lp_nm,
                                         DNASizeNanoMeters)

    NumAway = 2.5
    MinimumSeparationNanoMeters = NumAway*RadiusOfGyr
    # plots!
    fig = pPlotUtil.figure(figsize=(8,8))
    ax1 = plt.subplot(1,1,1)
    efficiency = 1./3
    plt.loglog(LoadConcentrationNgPerUl,MeanDistNano,'g-',linewidth=3,
               label="100% efficiency")
    # if we have lower efficiency, we are just lowering the number of molecules
    # available; singe the mean distance goes like 1/Sqrt(N), then if N_eff is
    # N/10 (meaning 10% efficiency), then the mean distance increases by
    # sqrt(10)~3
    LowerEfficiency= MeanDistNano*np.sqrt(1/efficiency)
    plt.loglog(LoadConcentrationNgPerUl,LowerEfficiency,
               'm.-',linewidth=3,
               label="{:d}% efficiency".format(int(efficiency*100)))
    plt.axvspan(xmin=MaxConcNgPerUl,xmax=plt.xlim()[-1],color='r',alpha=0.3,
                label="Impossible prep")
    plt.axhspan(ymin=plt.ylim()[0],ymax=MinimumSeparationNanoMeters,color='k',
                alpha=0.3)
    plt.axhspan(ymin=DNASizeNanoMeters,ymax=plt.ylim()[-1],color='k',alpha=0.3,
                label="Suboptimal for AFM")
    plt.fill_between(LoadConcentrationNgPerUl,y1=MeanDistNano,
                     y2=LowerEfficiency,
                     where=((MeanDistNano > MinimumSeparationNanoMeters) &
                            (LowerEfficiency < DNASizeNanoMeters)),
                     facecolor='k')
    plt.plot([], [], color='black', linewidth=15,label="Ideal Loading")
    # bit of a hack to get the label working
    pPlotUtil.lazyLabel("Concentration [ng/uL] (20uL Deposition)",
                        "Expected Mean DNA distance [nm]",\
                    "Deposit 20uL at a concentration of 0.2 ng/uL to 1ng/uL",
                        titley=1.1)
    # our circular DNA is roughly 607nm
    plt.axhline(DNASizeNanoMeters,linewidth=4,linestyle='--',
                label="L0={:.1f}".format(DNASizeNanoMeters))
    pPlotUtil.legend(loc='lower left',frameon=True)
    Molar= [NanoMolarLoaded[0],NanoMolarLoaded[-1]]
    pPlotUtil.secondAxis(ax1,"Molarity (nM)",limits=Molar,color="Blue",
                         secondY=False)
    plt.xlabel("Molarity (nM)")
    plt.tight_layout()
    fig.savefig("DepositionAdvice.png")
    """
    Figure 4/equation 7 from
    Latulippe, David R., and Andrew L. Zydney. 
    "Radius of Gyration of Plasmid DNA Isoforms from Static Light Scattering." 
    Biotechnology and Bioengineering 107, 2010
    """
    BasePairs = np.logspace(3.1,3.5)
    LengthsNm = BasePairs * NanometersPerBp
    # get the plot the DNA radii of gyration
    RadiusOfGyr = LinearRadiusOfGyration(Lp_nm,LengthsNm)
    RadiusOfGyrCirc = CircularRadiusOfGyration(Lp_nm,LengthsNm)
    NaiveRadiusOfGyr = np.sqrt(LengthsNm/Lp_nm) * Lp_nm / np.sqrt(6)
    # get the specific radii for our DNA for the model
    CircularSequenceRadiusCirc = CircularRadiusOfGyration(Lp_nm,
                                                          DNASizeNanoMeters)
    CircularSequenceRadiusLinear = LinearRadiusOfGyration(Lp_nm,
                                                          DNASizeNanoMeters)
    NumBasePairsCirc = KbpPerConstr*1000
    # plot everything
    StyleLin = dict(color='b',linewidth=2)
    StyleCirc = dict(color='r',linewidth=2)
    fig = pPlotUtil.figure()
    plt.plot(BasePairs,RadiusOfGyr,label="Linear WLC: {:d}nm".\
             format(int(CircularSequenceRadiusLinear)),**StyleLin)
    plt.plot(BasePairs,RadiusOfGyrCirc,
             label="Circular WLC: {:d}nm".\
             format(int(CircularSequenceRadiusCirc)),
             **StyleCirc)
    plt.axvline(NumBasePairsCirc,color='g',linewidth=2)
    XStart = plt.xlim()[0]
    ArrX = [XStart,NumBasePairsCirc]
    ArrY1 = [CircularSequenceRadiusCirc,CircularSequenceRadiusCirc]
    ArrY2 = [CircularSequenceRadiusLinear,CircularSequenceRadiusLinear]
    plt.plot(ArrX,ArrY1,linestyle="--",**StyleCirc)
    plt.plot(ArrX,ArrY2,linestyle="--",**StyleLin)
    pPlotUtil.lazyLabel("Base Pairs","Radius of Gyration (nm)",
                        "Circular DNA construct radius of gyration < 90nm",
                        frameon=True)
    pPlotUtil.savefig(fig,"Gyration.png")

if __name__ == "__main__":
    run()
