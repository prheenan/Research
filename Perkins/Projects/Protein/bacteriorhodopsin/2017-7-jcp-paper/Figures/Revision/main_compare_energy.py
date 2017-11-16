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

from Bio import SeqIO
from Bio.SeqUtils import seq1

def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    # get the
    record = SeqIO.read("B0R5N9.fasta", "fasta")
    seq_N_to_C = record.seq
    seq_C_to_N = seq_N_to_C[::-1]
    # get a key, start of the ED helix C->N pulling. The protein is listed N-C,
    # so we 
    key_C_to_N = "AKSTFGF"
    key_N_to_C = key_C_to_N[::-1]
    assert key_N_to_C in seq_N_to_C
    assert key_C_to_N in seq_C_to_N
    # POST: unique key
    idx_key = str(seq_C_to_N).index(key_C_to_N)
    seq_c_to_n = seq_C_to_N[idx_key:]
    # only look at the first N
    N_max = 160
    seq_c_to_n_measured = seq_c_to_n[:N_max]
    # determine the energies associated with transfer from water to lipid
    table = np.genfromtxt("./yamada_2016_energy.txt",defaultfmt="%s",dtype=str,
                          comments="#")
    table_str = "\n".join([",".join(t) for t in table])
    with open("table.csv","w") as f:
        f.write(table_str)
    three_letter_codes = table[:,0]
    one_letter_codes = [seq1(t) for t in three_letter_codes]
    energies_kJ_per_mol = [float(f) for f in table[:,1]]
    # make a dictionary key:value is aa code : energy of transfer
    aa_to_energy_dict = dict([ (k,v) for k,v in 
                               zip(one_letter_codes, energies_kJ_per_mol)])
    # determine the energy of the whole thing
    for s in aa_to_energy_dict:
        assert s in aa_to_energy_dict
    costs = [aa_to_energy_dict[s] for s in seq_c_to_n_measured]
    n_terminus_kJ_per_mol = 9.8
    total_cost_kJ_per_mol = sum(costs) + n_terminus_kJ_per_mol
    total_cost_kcal_per_mol = total_cost_kJ_per_mol/4.184
    cost_unfold_kcal_per_mol = -1 * total_cost_kcal_per_mol
    cost_unfold_kcal_per_mol_per_aa = cost_unfold_kcal_per_mol/N_max
    print("160 aa: {:s}".format(seq_c_to_n_measured))
    print("Total cost for {:d} aa: {:.1f} kcal/mol".\
          format(N_max,cost_unfold_kcal_per_mol))
    print("Total cost for per aa: {:.1f}".\
          format(cost_unfold_kcal_per_mol_per_aa))
                                                  

    


if __name__ == "__main__":
    run()
