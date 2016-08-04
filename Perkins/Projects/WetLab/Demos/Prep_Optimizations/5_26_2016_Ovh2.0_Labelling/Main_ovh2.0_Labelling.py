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

def LegendAndSave(Fig,SaveName,LegendLoc="upper right"):
    pPlotUtil.legend(loc=LegendLoc,frameon=True)
    pPlotUtil.savefig(Fig,SaveName,close=False)

def run():
    # get our actual sequence weight
    Primer = "AGAGTGGTCCTA"
    # get the sequence, tack on the overhang
    ProductStartLoc = 1606
    ProductEndLoc = 1794
    with open("mp13_plasmid_plasmid_seq.txt") as f:
        Seq = Primer + "".join([l for l in f.readlines()])
        Seq = Seq[ProductStartLoc:ProductEndLoc]
    # add in the overhang, note I also add an 'A' to represent the abasic site
    Overhang = "TAGGACCACTCT" + "A"
    Seq = Overhang + Seq
    SeqLen = len(Seq)
    MolecularWeight = Bio.SeqUtils.molecular_weight(Seq,
                                                    seq_type="DNA",
                                                    double_stranded=True,
                                                    circular=True)
    print("The Molecular weight of your {:d}bp sequence is {:.1f}kDa".\
          format(SeqLen,MolecularWeight/1000.))
    LengthDoubleStrandedBasePairs = len(Seq)
    Avogadro = 6.022e23
    MassGramsPerMolecule= MolecularWeight/Avogadro
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
    RadiusOfGyrLinearNanoMeters = LinearRadiusOfGyration(Lp_nm,
                                                         DNASizeNanoMeters)
    RadiusOfGyrNanoMeters = CircularRadiusOfGyration(Lp_nm,
                                           DNASizeNanoMeters)

    NumAway = 2.5
    MinimumSeparationNanoMeters = NumAway*RadiusOfGyrNanoMeters
    # plots!
    fig = pPlotUtil.figure(figsize=(8,8))
    ax1 = plt.subplot(1,1,1)
    efficiency = 1./5
    PerfectEfficiency = MeanDistNano
    plt.loglog(LoadConcentrationNgPerUl,PerfectEfficiency,'g-',linewidth=3,
               label="100% efficiency")
    # if we have lower efficiency, we are just lowering the number of molecules
    # available; singe the mean distance goes like 1/Sqrt(N), then if N_eff is
    # N/10 (meaning 10% efficiency), then the mean distance increases by
    # sqrt(10)~3
    LowerEfficiency= MeanDistNano*np.sqrt(1/efficiency)
    # N*L0 should be totally sufficient: DNA at site A must reach all the way
    # to site B (tether is complemtary region)
    LowerBoundDist = 1.5*DNASizeNanoMeters
    UpperBoundDist = 8*DNASizeNanoMeters
    plt.loglog(LoadConcentrationNgPerUl,LowerEfficiency,
               'm.-',linewidth=3,
               label="{:d}% efficiency".format(int(efficiency*100)))
    # bit of a hack to get the label working
    pPlotUtil.lazyLabel("Concentration [ng/uL] (20uL Deposition)",
                        "Expected Mean DNA distance [nm]",\
                        "Tuning Deposition Concentration to avoid dimerization",
                        titley=1.1)
    LegendAndSave(fig,"DepositionAdvice1.png")
    plt.axvspan(xmin=MaxConcNgPerUl,xmax=plt.xlim()[-1],color='r',alpha=0.3,
                label="Impossible prep")
    plt.axhspan(ymin=plt.ylim()[0],ymax=LowerBoundDist,color='k',
                alpha=0.3)
    plt.axhspan(ymin=UpperBoundDist,ymax=plt.ylim()[-1],color='k',alpha=0.3,
                label="Suboptimal for AFM")
    # our circular DNA is roughly 607nm
    plt.axhline(DNASizeNanoMeters,linewidth=4,linestyle='--',
                label="L0={:.1f}nm".format(DNASizeNanoMeters))
    LegendAndSave(fig,"DepositionAdvice2.png")
    IdealLoad = ((PerfectEfficiency > LowerBoundDist) &
                 (LowerEfficiency < UpperBoundDist))
    WhereIdealIdx = np.where(IdealLoad)
    plt.fill_between(LoadConcentrationNgPerUl,y1=LowerEfficiency,
                     y2=PerfectEfficiency,
                     where=IdealLoad,
                     facecolor='k')
    plt.plot([], [], color='black', linewidth=15,label="Ideal Loading")
    LegendAndSave(fig,"DepositionAdvice3.png")
    Molar= [NanoMolarLoaded[0],NanoMolarLoaded[-1]]
    pPlotUtil.secondAxis(ax1,"Molarity (nM)",limits=Molar,color="Blue",
                         secondY=False)
    plt.xlabel("Molarity (nM)")
    plt.tight_layout()
    pPlotUtil.savefig(fig,"DepositionAdvice4.png")
    """
    We load the 20uL into xuL of TE and heat; what volume should we use?
    """
    MaxMolarity = MolarityLoaded[WhereIdealIdx[0][0]]
    MinVolLog = np.log10(LoadVolumeMicroliters * 1e-6)
    MaxVolLog = np.log10(LoadVolumeMicroliters* 1e-6 * 100)
    DepositionVolumesLiters = np.logspace(MinVolLog,MaxVolLog)
    # get the total number of molecules, given we load XXuL
    TotalMolecules = (MaxMolarity*Avogadro)*(LoadVolumeMicroliters*1e-6)
    # get the spacing of these, as a function of the volume
    LoadVolumesMeters = DepositionVolumesLiters * 1e-3
    MeanVolume = LoadVolumesMeters/TotalMolecules
    MeanSpacingMeters = (MeanVolume * (3/(4*np.pi)))**(1/3)
    MeanSpacingNanoMeters = MeanSpacingMeters*1e9
    MeanSpacingRadiusOfGyr = MeanSpacingNanoMeters/RadiusOfGyrLinearNanoMeters
    # plot the mean spacing, in terms of radius of gyration, in terms of the
    # deposition volume
    fig = pPlotUtil.figure()
    NumRadii = 20
    ExtraMicroLiters = 30
    plt.plot(DepositionVolumesLiters * 1e6,MeanSpacingRadiusOfGyr,
             label="Mean spacing between DNA")
    plt.axvline(LoadVolumeMicroliters+ExtraMicroLiters,linewidth=3,color='r',
                label="20uL (Sample) + {:d}uL (TE) deposition".\
                format(ExtraMicroLiters))
    plt.axhline(NumRadii,linewidth=3,color='b',
                label="{:d} Radii of Gyration (Linear DNA)".format(NumRadii))
    plt.axhline(DNASizeNanoMeters/RadiusOfGyrLinearNanoMeters,
                linestyle="--",color='g',
                label="Contour Length")
    pPlotUtil.\
        lazyLabel("Microliters for Deposition",
                  "Mean Spacing (# of Linear Radii of Gyration)",
                  "With {:d}uL extra TE, DNA spaced >{:d} * Rg away".\
                  format(ExtraMicroLiters,NumRadii),
                  frameon=True)
    pPlotUtil.savefig(fig,"SpacingPerDeposition.png")
    """
    How long should we deposit for, assuming we want uniform coverage, taking
    into account the diffusion coefficient of the DNA?
    """
    """
    Diffusion Coefficient from Figure 1a (lower bound for 1um DNA), converted
    to m^2/s
    Robertson, Rae M., Stephan Laib, and Douglas E. Smith. 
    "Diffusion of Isolated DNA Molecules: Dependence on Length and Topology."
    PNAS, 2006
    """
    DiffusionCoeffMetersSquaredPerSec = 3*(1e-6)**2
    # 1mm to 100mm
    DistancesMeters = np.logspace(-5,-3)
    # use 2-D diffusion equation, worst-case
    TimeToDiffuse = lambda D,r : (r**2)/(4*D)
    TimesSeconds = TimeToDiffuse(DiffusionCoeffMetersSquaredPerSec,
                                 DistancesMeters)
    TimesHours = TimesSeconds/3600
    # write down the common surface depositions we use
    # 100uL deposition
    VolumeDepoLiters = 40e-6
    VolumeDepositions = [20e-6,40e-6]
    LabelDepositions = ["Mica","AFM"]
    # write down the depositions names, volumes (uL) and radii (m)
    # used, below, to plot
    Depositions = [["Mica",10e-6,9.5e-3/2],
                   ["AFM",80e-6,13e-3/2]]
    PlottingData = []
    for name,VolumeDepoLiters,RadiusMeters in Depositions:
        VolumeDepoMetersCubed = VolumeDepoLiters/1000
        # determine the effective average height of the bubble
        SurfaceAreaMetersSquared = np.pi * RadiusMeters**2
        HeightMeters = VolumeDepoMetersCubed/SurfaceAreaMetersSquared
        HeightMicrons = HeightMeters*1e6
        TimeDiffuseSecondsAFM = TimeToDiffuse(DiffusionCoeffMetersSquaredPerSec,
                                              HeightMeters)
        TimeDiffuseHoursAFM = TimeDiffuseSecondsAFM/3600
        Label = "{:s} Diffusion Time for {:d}uL deposition ~ {:.1f}hrs".\
                format(name,int(VolumeDepoLiters*1e6),TimeDiffuseHoursAFM)
        PlottingData.append([Label,HeightMicrons,TimeDiffuseHoursAFM])
    fig = pPlotUtil.figure()
    plt.plot(DistancesMeters*1e6,TimesHours)
    for Label,Height,Time in PlottingData:
        plt.plot(Height,Time,'ro',label=Label)
    pPlotUtil.lazyLabel("Distance (um)","Time to diffuse (hours)",\
        "DNA can't diffuse more than 0.5mm in a reasonable amount of time")
    pPlotUtil.savefig(fig,"DiffusionTimes.png")
    """
    Figure 4/equation 7 from
    Latulippe, David R., and Andrew L. Zydney. 
    "Radius of Gyration of Plasmid DNA Isoforms from Static Light Scattering." 
    Biotechnology and Bioengineering 107, 2010
    """
    BasePairs = np.logspace(2.3,3.5)
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
    ax = plt.subplot()
    plt.plot(BasePairs,RadiusOfGyr,label="Linear WLC",linestyle="--",**StyleLin)
    plt.plot(BasePairs,RadiusOfGyrCirc,label="Circular WLC",**StyleCirc)
    plt.axvline(NumBasePairsCirc,color='g',
                label="200bp, R={:d}nm".format(int(CircularSequenceRadiusCirc)),
                                             linewidth=2)
    LongerConstructBp = (3520-1607)+12+2
    LongerConstructNm = NanometersPerBp*LongerConstructBp
    LongerConstructRadiusNm = \
        int(CircularRadiusOfGyration(Lp_nm,LongerConstructNm))
    plt.axvline(LongerConstructBp,color='r',linestyle='-',
                label="{:d}bp, R={:d}nm".format(LongerConstructBp,
                                                LongerConstructRadiusNm),
                linewidth=2)
    SizeLimitsNanoMeters = \
        np.array([0,plt.xlim()[-1]])*NanometersPerBp
    XStart = plt.xlim()[0]
    ArrX = [XStart,NumBasePairsCirc]
    ArrY1 = [CircularSequenceRadiusCirc,CircularSequenceRadiusCirc]
    ArrY2 = [CircularSequenceRadiusLinear,CircularSequenceRadiusLinear]
    plt.plot(ArrX,ArrY1,linestyle="-.",**StyleCirc)
    plt.plot(ArrX,ArrY2,linestyle="-.",**StyleLin)
    pPlotUtil.lazyLabel("DNA Size, Base Pairs","Radius of Gyration (nm)",
    "Circular radius of gyration converges to linear for small constructs",
                        frameon=True,titley=1.1)
    pPlotUtil.secondAxis(ax,"Size (nm)",limits=SizeLimitsNanoMeters,
                         color="Blue",
                         secondY=False)
    pPlotUtil.xlabel("Size (nm)")
    pPlotUtil.savefig(fig,"Gyration.png")

if __name__ == "__main__":
    run()
