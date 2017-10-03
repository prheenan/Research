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

from Research.Perkins.Projects.Primers.Util import CommonPrimerUtil,IdtUtil
from GeneralUtil.python import PlotUtilities

def run():
    """
    Describes why we are picking the fluorophores we are
    """
    # get the excitation filter
    filter_file = "lumencor_631_28_nm.txt"
    arr = np.loadtxt(filter_file,skiprows=1)
    wavelength,filter_value = arr[:,0],arr[:,1]
    # convert from 0 to 1 transmission (from 0 to 100)
    filter_value /= 100
    # get the fluorophore
    # (see: http://www.fluorophores.tugraz.at/substance/418)
    fluorophore_file = "atto_5382.csv"
    arr = np.loadtxt(fluorophore_file,skiprows=1,delimiter=";",
                     usecols=(0,1,2,3))
    wavelength_fluorophore_excitation_nm,excitation = arr[:,0],arr[:,1]
    wavelength_fluorophore_emission_nm,emission = arr[:,2],arr[:,3]
    # get the emission filter 
    # see: laser2000.co.uk/semrock_filter.php?code=FF01-680/42-25
    emission_filter_file = "FF01-680_42.txt"
    arr_emit = np.loadtxt(emission_filter_file,skiprows=4)
    wavelength_emission_filter,emission_filter = arr_emit[:,0],arr_emit[:,1]
    # plot the fluorophores
    label_excite_emit = "\n({:s})".format(fluorophore_file)
    wavelength_limits_nm = [500,800]
    fig = PlotUtilities.figure((5,6))
    plt.subplot(2,1,1)
    plt.plot(wavelength,filter_value,
             label=("Filter (excitation)\n" +  filter_file))
    plt.plot(wavelength_fluorophore_excitation_nm,excitation,
             label="Fluorophore Excitation" + label_excite_emit,
             linestyle='--')
    PlotUtilities.lazyLabel("","Excitation efficiency","")
    plt.xlim(wavelength_limits_nm)
    PlotUtilities.no_x_label()
    plt.subplot(2,1,2)
    plt.plot(wavelength_emission_filter,emission_filter,
             label="Filter (emission) \n" + emission_filter_file)
    plt.plot(wavelength_fluorophore_emission_nm,emission,
             label="Flurophore Emission" + label_excite_emit,
             linestyle='--')
    plt.xlim(wavelength_limits_nm)
    PlotUtilities.lazyLabel("Wavelength (nm)","Emission efficiency","")
    PlotUtilities.savefig(fig,"./filter_comparisons.png")
    p1607F,_ = CommonPrimerUtil.Get1607FAnd3520R("../..")
    # add a dbco and a fluorophore...
    seq_full = [IdtUtil.Dbco5Prime()] + [s for s in p1607F] + \
               [IdtUtil.atto_633()]
    # get the IDT order
    opts = dict(Scale=IdtUtil.Scales._100NM,
                Purification=IdtUtil.Purifications.HPLC)
    order = IdtUtil.\
            SequencesAndNamesTuplesToOrder( [(seq_full,"AzideTestFluorophore")],
                                            **opts)
    IdtUtil.PrintAndSave(order,"./test_azide_fluorophore.txt")
    

if __name__ == "__main__":
    run()
