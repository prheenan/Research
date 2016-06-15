from collections import Counter
# product gives all possible...
from itertools import product
import numpy as np
DEF_BASES_DNA = ["A","C","G","T"]

def SanitizeSequence(Seq):
    """
    Sanitize a sequence string so that all the methods can use it 

    Args:
        Seq: The sequence to santize
    Returns:
        The sanitized Sequence
    """
    return Seq.strip().upper()

def Complement(strV):
    """
    Gets the complement, dna-wise of strV

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        strV: the complement of the DNA to get
        lenUniKmer: the length of the kmer to parse
    
    Returns:
        the complement of the sequence
    """    
    dictDNA = {'A':'T',
               'T':'A',
               'G':'C',
               'C':'G'}
    return ''.join([ dictDNA[base] for base in strV ])

def ReverseComplement(strV):
    """
    Gets the reverse complement of a string, DNA-wise

    Args:
        strV: the reverse complement of the DNA to get

    Unit Tested by MainTestUtil.TestKmerUtil    

    Returns:
        the reverse complement of the sequence
    """   
    return Complement(strV[::-1])

def Kmers(rawStr,lenUniKmer):
    """
    Gets the in-order occurences (*not* a set)
    of all kmers of length lenUniKmer in rawStr

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        rawStr: The string to parse
        lenUniKmer: the length of the kmer to parse
    
    Returns:
        the list of all kmers in the string
    """        
    plasmidLen = len(rawStr)
    # get all the kmers
    uniqueKmers = [rawStr[i:i+lenUniKmer] 
                   for i in range(plasmidLen-lenUniKmer+1)]
    return uniqueKmers

def KmerSet(rawStr,lenUniKmer):
    """
    Gets the *set* of all kmers in rawStr

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        rawStr: The string to parse
        lenUniKmer: the length of the kmer to parse
    
    Returns:
        the set of all kmers in the string
    """     
    return set(Kmers(rawStr,lenUniKmer))

def GetWraparoundAndReverseComplement(rawStr,lenK):
    """
    Given a raw string and a kmer length, gets the string with 'wrapraround'
    and reverse complements (also with wrapraround). Wraparound is due to
    (e.g.) a circular piece of DNA

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        rawStr: the raw string to get the wraparound/revcomp of
        lenK: the length of the kmers (== length of wraparond)
    Returns:
        tuple of fullForward/fullbackwards string

    """
    # function to add 'wraparound' to a string
    AddWraparound = lambda StrV,lenV : StrV[-lenV:] + StrV + StrV[:lenV]
    # get the full foward string
    FullForwardStr = AddWraparound(rawStr,lenK)
    # get the Full backwards string, just the reverse complement of the forward
    FullBackwardStr = ReverseComplement(FullForwardStr)
    return FullForwardStr,FullBackwardStr

def CircularizedKmers(rawStr,lenKmers):
    """
    Gets all the kmers in the forward and reverse sequences, accounting 
    for a circular plasmid (ie: accounting for moving across the two ends

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        rawStr: The raw string to use
        lenKmers: the unique kmer size
    Returns:
        set of kmers in the forward, reverse srings, including overlap
    """
    FullForwardStr,FullBackwardStr = GetWraparoundAndReverseComplement(rawStr,
                                                                       lenKmers)
    # return all the kmers in either set
    return KmerSet(FullForwardStr,lenKmers) | \
        KmerSet(FullBackwardStr,lenKmers)

def AllPossibleKmers(lenV,bases=DEF_BASES_DNA):
    """
    Gets all Possible Kmers of length lenV. Only recommended for small lengths

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        lenV: length of the string to get
        bases: possible bases
    Returns:
        all possible kmers of the given length.
    """
    # cartesian product gives all possible combinations
    return set([''.join(s) for s in product(bases,repeat=lenV)])
        
def GetKmersNotAppearing(plasmidStr,lenUniKmer,bases=DEF_BASES_DNA):
    """
    Gets all kmers of length 'lenUniKmer' not appearing in plasmidStr or its
    reverse complement, also accounting for 'wraparound' on the plasmid.

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        plasmidStr: string to search
        lenUniKmer: length of the kmers to search for
    Returns:
        set of all kmers of the given length *not* appearing in the input string
    """
    # print: plasmid is just a sequence, no whitespace
    # dont want any overlaps in the string, its reverse, its reverse complement,
    # or its complement
    uniqueKmers = CircularizedKmers(plasmidStr,lenUniKmer)
    # POST: have the sequence for every unique k-mer in the plasmid.
    # remove the set-wise intersection of all
    # possible k-mers and the unique ones in the plasmid. 
    allPossible = AllPossibleKmers(lenUniKmer,bases=bases)
    # This gives all the kmers which are 'left' (ie: non overlapping)
    nonOverlapping = allPossible - uniqueKmers
    return nonOverlapping

def GetKmerCounts(InputStr,kmerLen):
    """
    Counts the occurences of all the kmers in inputStr, returns as a dictionary

    Unit Tested by MainTestUtil.TestKmerUtil

    Args:
        InputStr: String to search
        kmerLen: length of kmers to search for
    Returns:
        in-order count for each kmer to count
    """
    # first, get all the kmers
    occuringKmers = Kmers(InputStr,kmerLen)
    # make a counter for the kmers
    return Counter(occuringKmers)


def BestKmer(plasmidSeq,k):
    """
    Get the best kmer to add, assuming all kmers appear: which is the least
    likely? Looks for overhangs and in the reverse complement

    Args:  
        plasmidSeq: the sequence to look for
        k: the k to use
    Returns:
        tuple of <kmer,occurences>
    """
    # get *least* common item occuring at least once
    # (see
    # https://docs.python.org/3.2/library/collections.html#collections.Counter)
    return LeastOccuringKmers(plasmidSeq,k,N=1)

def LeastOccuringKmers(plasmidSeq,k,N=None):
    """
    Get the N least occuring kmers are <kmer:value> pairs
    Args:  
        plasmidSeq: the sequence to look for
        k: the k to use
        N: the kmers to use. If none, 
    Returns:
        list of tuples of <kmer,occurences>, ordered by least occurening first
    """
    mostCommon = AllKmerCounts(plasmidSeq,k).most_common()[::-1]
    if (N is None):
        # return them all
        return mostCommon
    else:
        # return just the N we want
        return mostCommon[:N]
        
def AllKmerCounts(plasmidSeq,k):
    """
    Looks for all kmers in the forward, reverse complement, and overhanging
    regions of length k
    
    Args:
        plasmidSeq: the sequence to look for
        k: the k to use
    Returns:
        Python Counter (dictionary) object with key:values of <kmers:counts>
    """
    FullForwardStr,FullBackwardStr = \
            GetWraparoundAndReverseComplement(plasmidSeq,k)
    counterFwd= GetKmerCounts(FullForwardStr,k)
    counterBack= GetKmerCounts(FullBackwardStr,k)
    # get the total count of the kmers occuring 
    totalCount = counterFwd + counterBack
    return totalCount


def AssertValidKmers(plasmidStr,kmerList):
    """
    Make sure none of the kmers on the list are present in plasmidStr 

    Unit Tested by 'TestKmers.TestKmerPrimersValid'

    Args:
        plasmidStr: the plasmid to read in
        kmer: list of kmers to look through. should *not* be in plasmidStr
    """
    # get all the kmers associated with the plasmid
    lenV = len(kmerList)
    # get the length of a single kmer. probably ineffieicent
    kmerLen = len(list(kmerList)[0])
    # get all the kmers in the plasmid
    kmersAppearingInPlasmid = CircularizedKmers(plasmidStr,kmerLen)
    intersection = (kmersAppearingInPlasmid & set(kmerList))
    # test that the kmers indeed arent in the plasmid
    assert len(intersection) == 0, "{:s}".format(intersection)


def AssertPrimersValid(plasmidStr,kmerLen,ListOfPrimers):
    """
    Given a primer, plasmidStr, and kmer, ensures that none of the kmers
    in the primer are in the plasmid. Useful for something like an overhang.

    Unit Tested by 'TestKmers.TestKmerPrimersValid'

    Args:
        plasmidStr: plasmid string to search in
        kmerLen: length of kmers we dont want
        listOfPrimers: list of all the primers we want to check
    """
    primerKmers = set([k for p in ListOfPrimers
                       for k in Kmers(p,kmerLen)])
    AssertValidKmers(plasmidStr,primerKmers)

def SmallestNonCoveredKmer(primerStr,thresh=1,CheckCircular=False):
    """
    Gets the smallest k (k*) such that at least 'thresh' (k*)mers are not
    present in primer string, or its reverse complement, or circular region

    This is useful if we want to find a k* for overhangs (ie: kmers of length
    k* will be the smallest 'alphabet' we should use, guarenteeing no overlap
    of length (k*-1)

    Args:
        primerStr: primer string to check
        thresh: need to have at least this many kmers unaccounted for. 
        CheckCircular: if true, checks the sequence, its reverse complement,
        and its 'wraparound' (ie: assumes a plasmid)
    """
    n = len(primerStr)
    numBases = 4
    minK = int(np.ceil(np.log(thresh)/np.log(numBases)))
    for k in range(minK,n):
        if (CheckCircular):
            unique= CircularizedKmers(primerStr,k)
        else:
            unique = set(Kmers(primerStr,k))
        numPossible = (numBases)**k
        # check if we found the correct k
        if ((numPossible-len(unique)) >= thresh):
            break
    return k
