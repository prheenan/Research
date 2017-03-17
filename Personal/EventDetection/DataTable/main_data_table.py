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
sys.path.append("../../../../")

from Research.Personal.EventDetection.Util import Plotting,InputOutput
from GeneralUtil.python import CheckpointUtilities,GenUtilities,PlotUtilities

class data_table:
    def __init__(self,category):
        self.sizes =[d.Force.size for d in category.data]
        self.n_events = [len(d.Events)for d in category.data]
        self.meta_dicts = [d.Meta for d in category.data]
        self.velocity = category.velocity_nm_s
        self.directory = category.directory 
    def n_fec(self):
        return len(self.n_events)
    def n_total_events(self):
        return sum(self.n_events)
    def counts(self,list_of_counts):
        to_ret = []
        for num in list_of_counts:
            to_ret.append(self.n_events.count(num))
        return to_ret
    def counts_ge(self,num):
        to_ret = []
        idx_ge = (np.where(np.array(self.n_events) >= num))[0]
        return idx_ge.size
    def get_relevant_numbers(self,count_numbers,count_ge):
        zero_idx = np.where(np.array(self.n_events) < 0.5)[0]
        assert zero_idx.size == 0 , "Missing an event"
        n_small = self.counts(count_numbers)
        n_ge = self.counts_ge(count_ge)
        assert (sum(n_small) + n_ge) == self.n_fec()
        mean_n_points = np.mean(self.sizes)
        std_n_points = np.std(self.sizes)
        return [self.velocity,self.n_fec(),mean_n_points,std_n_points] + \
            n_small + [n_ge]


def get_data_table(cache_directory,limit,force_read):
    positives_directory=  InputOutput.get_positives_directory()
    categories = InputOutput.get_categories(positives_directory)
    force_read=False
    categories = InputOutput.read_categories(categories,force_read,
                                             cache_directory,limit)
    table_rows = []
    for c in categories:
        table_rows.append(data_table(c))
    return table_rows
        
    
def run():
    """
    <Description>

    Args:
        param1: This is the first param.
    
    Returns:
        This is a description of what is returned.
    """
    cache_directory = "../_1ReadDataToCache/cache/"
    limit = 200
    force_read=False
    table_rows= \
        CheckpointUtilities.getCheckpoint(cache_directory + "table.pkl",
                                          get_data_table,force_read,
                                          cache_directory,limit,
                                          force_read)
    count_numbers = [1,2,3]
    count_ge = [4]         
    file_out = cache_directory + "table.txt"
    formatted_string = ""
    round_integer = lambda x : int(np.round(x))
    round_decimal = lambda x,d : int(np.round(x,d))
    round_thousands = lambda x: round_decimal(x,-3)
    rounding_funcs = [round_integer,round_integer,round_thousands,
                      round_thousands,round_integer,round_integer,round_integer,
                      round_integer]
    numbers = []
    # get the number for the data table
    for t in table_rows[::-1]:    
        numbers.append(t.get_relevant_numbers(count_numbers,count_ge))
    # round them
    formatted_numbers = []
    for nums in numbers:
        fmt_tmp = ["{:d}".format(f(n)) for f,n in zip(rounding_funcs,nums)]
        formatted_numbers.append(fmt_tmp)
    # make a latex table 
    join_str = " & "
    event_headers = [(r"N$_{\mathrm{e}=" + (" {:d}").format(n) + r"}$" )
                     for n in count_numbers]
    event_headers += [(r"N$_{\mathrm{e}\ge"+("{:d}").format(count_ge[0]))+"}$"]
    headers = ["v [nm/s]",r"N$_\mathrm{curves}$","$\mu_{\mathrm{Curve Size}}$",
                "$\sigma_{\mathrm{Curve Size}}$"] + event_headers
    end_str = r"  \\ \hline" + "\n"
    format_str = join_str.join(headers) + end_str
    for f in formatted_numbers:
        format_str += (join_str.join(f) + end_str)
    with open(file_out,'w') as f:
        f.write(format_str)

if __name__ == "__main__":
    run()
