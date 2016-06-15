# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys
import networkx


def GetPathThroughKmers(kmerGraph):
    """
    Given a KmerGraph (see "GetKmerGraph"), finds the longest cycle in the 
    graph, in an attempt to find the longest path. 

    Args:  
        kmerGraph: see GetKmerGraph
    Returns:
        longest list of kmers which have no 'bad' overlaps (according to edges)
    """
    # get all the simple paths between every node and every other node
    nodes = kmerGraph.nodes()
    allCycles = list(networkx.simple_cycles(kmerGraph))
    # sort the cycles by their length
    cycleLengths = [len(c) for c in allCycles]
    possibleCycles = [allCycles[idx] for idx in np.argsort(cycleLengths)][::-1]
    # start off with the maximum cycle...
    fullCycle = possibleCycles[0]
    # add in the last element to make it truly a cycle.
    fullCycle.append(fullCycle[0])
    return fullCycle

def GetKmerGraph(validAlphabet):
    """
    Given an 'alphabet' of subtrings (e.g. kmers), returns a graph where
    each node is a substring, and each edge represents 'the overlap between
    these two substrings is in the alphabet'

    Args:  
        validAlphabet: list of valid substrings
    Returns:
        XXX
    """
    # first, construct a graph. Each node is a substring, each edge represents
    # 'the overlap between these is also a substring'
    graph = networkx.DiGraph()
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
