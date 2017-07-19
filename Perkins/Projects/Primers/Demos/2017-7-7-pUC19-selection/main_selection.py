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

from Bio import SeqIO,Restriction 

def run():
    """
    """
    file_name = "addgene-plasmid-50005-sequence-74677.gbk"
    seqs = [s for s in SeqIO.parse(file_name, "genbank")]
    seq_record = seqs[0]
    # get the sequence as a string
    seq_as_string = str(seq_record.seq[:])
    # write down the different restriction enzymes
    re_sequences = [ Restriction.EcoRI,
                     Restriction.SacI,
                     Restriction.KpnI,
                     Restriction.SmaI,
                     Restriction.BamHI,
                     Restriction.XbaI,
                     Restriction.SalI,
                     Restriction.SbfI,
                     Restriction.PstI,
                     Restriction.SphI,
                     Restriction.HindIII]
    for re_obj in re_sequences:
        seq = re_obj.site
        name = str(re_obj)
        full_circular = seq_as_string + seq_as_string[:int(len(seq)/2)+1]
        count = full_circular.count(seq)
        print("{:s} ({:s}) has {:d} occurences".format(name,seq,count))
        assert count == 1 , "Restriction enzyme found not exactly once"
    # POST: all restriction enzymes found exactly once
    # determine:
    # (0) does it digest plasmids well?
    # https://www.neb.com/~/media/NebUs/Files/Brochures/RestEndo_TechGuide.pdf
    # 'Time-Saver Qualified'
    # (1) can be heat inactivated? 
    # www.neb.com/tools-and-resources/usage-guidelines/heat-inactivation
    
    


if __name__ == "__main__":
    run()
