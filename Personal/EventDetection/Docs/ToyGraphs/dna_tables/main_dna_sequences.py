# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

def read_seq(input_file):
    with open(input_file) as f:
        # first line is comments
        f.readline()
        # second line is sequence
        seq = f.readline()
    return seq

def format_seq(seq,
               terminal_biotin="B",
               terminal_DBCO="D",
               internal_biotin="B"):
    seq = seq.replace("5dbcoteg",terminal_DBCO)
    seq = seq.replace("3bioteg",terminal_biotin)
    seq = seq.replace("dbco",terminal_DBCO)
    seq = seq.replace("5biotin",terminal_biotin)
    seq = seq.replace("ibiodt",internal_biotin)
    seq = seq.replace("/","")
    return seq

class seqs:
    def __init__(self,fwd_primer,rev_primer,hairpin):
        sanit = lambda x: x.replace(" ","").lower()
        self.fwd = sanit(fwd_primer)
        self.rev = sanit(rev_primer)
        self.hairpin = sanit(hairpin)
    def fwd_rev_hairpin_formatted(self):
        return format_seq(self.fwd),\
            format_seq(self.rev),\
            format_seq(self.hairpin)

def get_sequences(base_dir):
    return seqs(read_seq(base_dir  + "1607F_DBCO"),
                read_seq(base_dir  + "3520R_4xBio"),
                read_seq(base_dir  + "68nt_hairpin"))

def get_latex_table(sequences):
    formatted_seqs = sequences.fwd_rev_hairpin_formatted()
    line_end = " \\\\ \e \n"
    to_ret = "Name & Sequence" + line_end
    names = ["Forward primer for 650nm DNA",
             "Reverse 650nm DNA",
             "68nt hairpin"]
    rows = ["{:s} & {:s}".format(name,seq) 
            for name,seq in zip(names,formatted_seqs)]
    to_ret += line_end.join(rows) +  line_end
    return to_ret


def run(base_dir="./"):
    """
    """
    sequences = get_sequences(base_dir=base_dir)
    table = get_latex_table(sequences)
    print(table)
    
if __name__ == "__main__":
    run()
