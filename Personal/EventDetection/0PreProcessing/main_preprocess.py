# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys,os

sys.path.append("../../../../")
import csv,re

from Research.Perkins.AnalysisUtil.ForceExtensionAnalysis import FEC_Util
from GeneralUtil.python import GenUtilities,CheckpointUtilities
from IgorUtil.PythonAdapter import TimeSepForceObj,DataObj,ProcessSingleWave



def get_events(file_name):
    """
    given a file formatted like an event file, reads them in 

    Args:
        file_name: path to read.
    Returns:
        list of tupels like <name of event, start index, end of index>
    """
    # get the file without the '.pxp' extension
    f_no_ext = file_name[:-4]
    # add on the suffix for the events
    f_events = f_no_ext + "_events.txt"
    # read in the (csv) file
    assert os.path.isfile(f_events) , \
        "Couldn't find event file for {:s}".format(f_events)
    # POST: event file exists
    # the syntax is <name,start of event, end of event> separated by commss
    with open(f_events) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        data = [ [str(r[0].strip()),int(r[1]),int(r[2])]
                for r in spamreader]
    return data

def read_single_directory_with_events(directory):
    """
    Reads in pxp files and their associated events

    Args:
        directory: where to read the pxp (each is assumed to have 
    """
    pxp_files,data = FEC_Util.read_single_directory(directory)
    # each pxp file should have events associated with it 
    events = [get_events(f) for f in pxp_files]
    # create a tuple of everything...
    combined = [(file_name,dat,ev)  \
                for file_name,dat,ev in zip(pxp_files,data,events)]
    return combined

def get_id(x):
    """
    Gets the id associated with the string x

    Args:
        x: probably something like the name of a trace
    Returns:
        id (e.g. Image0004_Force would give Image0004)
    """
    id_regexpr = r"""
                 (\D+ # any number of non-numbers
                 \d+) # any number of digits
                 \D"""  # anything *not* a digit
    match = re.match(id_regexpr,x,re.VERBOSE)
    assert match is not None , "Couldn't find id of {:s}".format(x)
    return match.group(1)

def set_events_of_data(data,events):
    """
    sets the events of all the data traces

    Args:
        data: list of time_sep_force objects
        events: list of <Trace Name, start,end> giving the events
    Returns:
        nothing but sets the object's events appropriately
    """
    id_data = [get_id(d.Meta.Name) for d in data]
    id_events = [get_id(e[0]) for e in events]
    # determine matches; may have multiple events
    eq = lambda x,y: x.strip().lower() == y.strip().lower()
    for id_data_tmp,d in zip(id_data,data):
        # find which index (in id_events) corresponds to id_data_tmp
        matching_idx = [j for j,id_ev in enumerate(id_events) 
                        if eq(id_ev,id_data_tmp)]
        # get the actual events
        events_matching = [events[i] for i in matching_idx]
        # add the events to the TimeSepForce Object. Note that
        # an event 'e' is like <name,start,end> so we just get the starts 
        # and ends
        starts_and_ends = [e[1:] for e in events_matching]
        Events = [TimeSepForceObj.Event(*s) for s in starts_and_ends]
        # set the events of the data
        d.set_events(Events)
    # POST: all events have been set

def run():
    """
    
    """
    base_directory="/Volumes/group/4Patrick/CuratedData/"
    output_base_directory = base_directory + "Masters_CSCI/"
    # copy each of the input directories into CSV files in the output directory
    dna_relative = "CircularDNA/ovh2.0-1p9kbp/ForceExtensionCurves/RawPxpFiles/"
    relative_input_dir = [dna_relative + "Scratch/"]
    absolute_input_dir = [(base_directory + d) for d in relative_input_dir]
    files_data_events = [read_single_directory_with_events(d) 
                         for d in absolute_input_dir]
    for i,d in enumerate(absolute_input_dir):
        # go through each PXP in this directory...
        for file_name,data,ev in files_data_events[i]:
            set_events_of_data(data,ev)
            # POST: all data are set. go ahead and save them out.
            for dat in data:
                output_path = d + "tmp.csv"
                FEC_Util.save_time_sep_force_as_csv(output_path,dat)
                time_sep_force = FEC_Util.\
                    read_time_sep_force_from_csv(output_path,has_events=True)
    

if __name__ == "__main__":
    run()
