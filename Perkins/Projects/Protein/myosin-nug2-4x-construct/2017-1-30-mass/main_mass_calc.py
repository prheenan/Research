# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy import constants

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    # note: our construct is 4x NUG2
    # see: http://www.rcsb.org/pdb/explore.do?structureId=1mi0
    molecular_weight_daltons = 14712.24 * 4
    dalton_in_kg = constants.physical_constants['atomic mass constant'][0]
    avogadro_atoms_per_mol = constants.\
        physical_constants['Avogadro constant'][0]
    molecular_weight_kg= molecular_weight_daltons * dalton_in_kg
    molecular_weight_ug = molecular_weight_kg * 1e3 * 1e6
    # write down the concentration we know
    concentration_ug_per_mL = 0.4/200 * 1e3
    conventration_molecules_per_mL = \
        concentration_ug_per_mL/molecular_weight_ug
    # convert to moles per mL
    concentration_moles_per_mL = \
        conventration_molecules_per_mL/avogadro_atoms_per_mol
    # convert to molarity (moles per L), 1000mL/L
    concentration_molarity = concentration_moles_per_mL * 1e3
    # convert to nanomolar
    concentration_nM = concentration_molarity*1e9
    out = "The deposited concentration of the 4x NUG construct is: {:d}nM".\
          format(int(np.round(concentration_nM)))
    print(out)

if __name__ == "__main__":
    run()
