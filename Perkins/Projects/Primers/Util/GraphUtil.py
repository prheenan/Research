# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import networkx as nx
import KmerUtil

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
    PathLength = PathLengthMin + 1
    PrimerSet = set()
    for i in vertices:
        for j in vertices:
            RawPaths = nx.all_simple_paths(Graph, source=i,
                                           target=j, cutoff=PathLength)
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


def GetKmerGraph(validAlphabet):
    """
    Given an 'alphabet' of subtrings (e.g. kmers), returns a graph where
    each node is a substring, and each edge represents 'the overlap between
    these two substrings is in the alphabet'

    Args:  
        validAlphabet: list of valid substrings
    Returns:
        NetworkX.DiGraph (directed graph). Each node is a valid symbol from
        the alphabet of kmers, each edge represents a 'safe' path, all over
        lapping kmers are in the alphabet (ie: are 'safe')
    """
    # first, construct a graph. Each node is a substring, each edge represents
    # 'the overlap between these is also a substring'
    graph = nx.DiGraph()
    alphaSet = set(validAlphabet)
    for u in alphaSet:
        for v in alphaSet:
            lenU = len(u)
            lenV = len(v)
            minLen = min(lenU,lenV)
            works=True
            # loop through, determine if the overlap between these kmers
            # is in the alphabet. 
            for i in range(1,minLen):
                if (u[i:] + v[:i]) not in validAlphabet:
                    works=False
                    break
            if (works):
                graph.add_edge(u,v)
    return graph
