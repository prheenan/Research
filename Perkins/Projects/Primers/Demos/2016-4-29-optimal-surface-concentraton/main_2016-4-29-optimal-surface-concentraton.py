# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys



def run():
    """
    Given the dsDNA we want to produce, calculates the optimal concentration
    for circularization
    """
    # amplify 3520R+1607F
    length_bp = (3520-1607)
    # get the volume associated with this length, assuming .34nm/bp
    bp_to_nm = 0.34
    length_nm = bp_to_nm * length_bp
    length_m = length_nm * 1e-9
    # get the radius of gyration, a=(one persistence length), from
    # https://en.wikipedia.org/wiki/Radius_of_gyration
    a = 43e-9
    N = int(np.ceil(length_m/a))
    length_gyrated_m = (1/np.sqrt(6)) * np.sqrt(N) * (a)
    # assume DNA is a sphere, what volume does it take up?
    vol_m_cubed = (4*np.pi/3)*(length_gyrated_m)**3
    # what concentration is this?
    critical_conc_molecules_per_m_cubed = 1./vol_m_cubed
    Avo = 6.02*1e23
    # 1m^3/1000L
    critical_conc_molecules_per_L = critical_conc_molecules_per_m_cubed *\
                                    1/1000
    # write the molarity, for kicks (moles per L)
    critical_molarity_moles_per_L =  critical_conc_molecules_per_L/Avo
    # 1L/(10^(6) mu L)
    critical_conc_molecules_per_uL = critical_conc_molecules_per_L*\
                                     1/1e6
    """
    Now we need to figure out what ng/uL that concentration corresponds to

    For data in DNA weight and concntrations, see (aslso saved as pdf here)
    www.neb.com/tools-and-resources/usage-guidelines/nucleic-acid-data
    """
    # mass of dna basepair in kg, converted from daltons
    dalton_per_bp = 650
    kg_per_dalton = 1.66*1e-27
    kg_per_bp = kg_per_dalton*dalton_per_bp
    # get the mass for all the DNA
    kg_per_molecule = length_bp * kg_per_bp
    # get the critical concentration in kg/m^3
    conc_kg_per_uL_cubed = critical_conc_molecules_per_uL * \
                           kg_per_molecule
    # 1ng/(10^(-12) kg)
    conc_ng_per_uL_cubed =  conc_kg_per_uL_cubed/(1e-12)
    # roughly speaking by NEB, 1ng/uL ~ 0.15mM/50 = 0.003mM
    print(("Gyrated Length:{:.2g}nm\n" +
          ("Critical Conc :{:.2g} molecules/uL = {:.3g}mM\n")+
          ("Critical Conc :{:.2f} ng/uL\n")).\
          format(length_gyrated_m*1e9,
                 critical_conc_molecules_per_uL,
                 1000*critical_molarity_moles_per_L,
                 conc_ng_per_uL_cubed))
    

if __name__ == "__main__":
    run()
