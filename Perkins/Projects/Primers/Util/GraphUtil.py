# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import networkx as nx
import KmerUtil

def FilterGraphByEdges(Graph,FilterFunction):
    """
    Given a graph and a function, creates a new graph by applying a boolean
    filter to each edge
    
    Args:
        Graph: networkx.graph to use
        FilterFunction: function taking in (i,j), where is one edge, j is 
        another. In DiGraph, directed from i to j
    Returns:
        filtered version of graph
    """
    ToRet = Graph.copy()
    ToRet.clear()
    for (u,v) in Graph.edges_iter():
        if FilterFunction(u,v):
            ToRet.add_edge(u,v)
    return ToRet
            
def GetValidPrimersFromGraph(Graph,PrimerLength):
    """
    Given a kmer graph (see GetKmerGraph), returns all valid primers of length
    DesiredPrimerLength, using edges as links between valid kmers

    Args:
        Graph: output of GetKmerGraph
        PrimerLength: How long the ending primers must be
    Returns: 
        list of unique, valid primers  of length DesiredPrimerLength
    """
    edges = Graph.edges()
    KmerLen = len(edges[0][0])
    vertices = Graph.nodes()
    # determine the maximum number of edges we need
    NumBpPerNode = KmerLen
    # first node is 'free'; how many more after that do we need?
    # this is the (P-1), where P is the minimum length of the path
    NumEdgesNeeded = np.ceil((PrimerLength-(NumBpPerNode*1))/NumBpPerNode)
    PathLengthMin = NumEdgesNeeded +  1
    # we add one to PathLengthMin for the cutoff, to avoid edge effects
    PathLength = int(PathLengthMin + 1)
    PrimerSet = set()
    for i,u in enumerate(vertices):
        print(i)
        for j,v in enumerate(vertices):
            RawPaths = list(nx.all_simple_paths(Graph, source=u,
                                                target=v, cutoff=PathLength))
            # Get the paths of the appropriate length
            LongEnoughPaths = [p for p in RawPaths if len(p) >= PathLengthMin]
            # convert the paths to strings; each string is a valid 'psuedo
            # plasmid'
            PsuedoPlasmid = ["".join(p) for p in LongEnoughPaths]
            # get the valid primers in the plasmid by breaking into
            # PrimerLen-mers
            ValidPrimers = [primer for p in PsuedoPlasmid
                            for primer in KmerUtil.Kmers(p,PrimerLength)]
            # add to the set of primers
            NewSetOfPrimers = set(ValidPrimers)
            PrimerSet = PrimerSet | NewSetOfPrimers
    # POST: all done, return the set of primers as list
    return sorted(list(PrimerSet))


def GetKmerGraph(validAlphabet,ValidFunc = lambda x : True):
    """
    Given an 'alphabet' of subtrings (e.g. kmers), returns a graph where
    each node is a substring, and each edge represents 'the overlap between
    these two substrings is in the alphabet'

    Args:  
        validAlphabet: list of valid substrings
        ValidFunc : by default, if every substring is in the alphabet, we 
        add. If non-default, we apply this function to the two edges, as well.
        It takes in a single primer, returns true if it is valid
    Returns:
        NetworkX.DiGraph (directed graph). Each node is a valid symbol from
        the alphabet of kmers, each edge represents a 'safe' path, all over
        lapping kmers are in the alphabet (ie: are 'safe')
    """
    # first, construct a graph. Each node is a substring, each edge represents
    # 'the overlap between these is also a substring'
    graph = nx.DiGraph()
    alphaSet = set(validAlphabet)
    LenSingle = len(list(validAlphabet)[0])
    for i,u in enumerate(alphaSet):
        for j,v in enumerate(alphaSet):
            works=True
            # loop through, determine if the overlap between these kmers
            # is in the alphabet.
            NewSet = KmerUtil.KmerSet(u+v,LenSingle)
            AllValid = ValidFunc(str(u) + str(v))
            if (len(NewSet - alphaSet) == 0) and AllValid:
                graph.add_edge(u,v)
    return graph
