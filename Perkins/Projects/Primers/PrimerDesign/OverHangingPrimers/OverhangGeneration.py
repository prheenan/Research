# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

import Util.KmerUtil as KmerUtil
import Util.EmbossUtil as EmbossUtil
import Util.MeltingTemperatureUtil as MeltingUtil
from Util.EmbossUtil import Primer
import Util.GraphUtil as GraphUtil
import Util.AlignUtil as AlignUtil
import Util.IdtUtil as IdtUtil
import Util.python.CheckpointUtilities as CheckUtil
import networkx
import Util.CommonPrimerUtil as CommonPrimerUtil

IDT_MIN_NUM = 15


class OverhangPrimerInfo:
    """
    class used to keep track of all the information we get when we try to 
    generate primers
    """
    def __init__(self,Plasmid,PossibleKmers,KmerGraph,AllPrimers):
        """
        Args:
            Plasmid: the raw string of the plasmid

            Possible Kmers: all possible kmers which we can choose from
            (ie: which don't overlap with the plasmid or its reverse complement
    
            KmerGraph: nodes are kmers, edges represent a 'safe' path between
            two kmers (ie: the overlap between the kmers isn't found in the 
            plasmid
         
            AllPrimers: all the possible primers we obtained
            
        """
        self.Plasmid = Plasmid
        self.PossibleKmers = PossibleKmers
        self.KmerGraph = KmerGraph
        self.AllPrimers = AllPrimers
        # set later...
        self.BestPrimers = None
        # alphabetK is the size of the k-mers we use to construct the plasmid
        self.AlphabetK = len(list(PossibleKmers)[0])

    def SetBestPrimers(self,BestPrimers):
        """
        Sets the best primers (presummably, not self-dimerizing with a good
        melting temperature) to BestPrimers:
        Args:
            BestPrimers: list of primer objects
        """
        self.BestPrimers=BestPrimers
    def SetAllPossiblePrimers(self,PossiblePrimers):
        """
        Sets the instance variable AllPrimers to the possible primers, 
        by virtue of having no k-mers in common with the plasmid
        
        Args:
            PossiblePrimers: list of all primers
        """
        self.AllPrimers = PossiblePrimers
    
def GetValidKmers(plasmidStr,kmer):
    """
    Gets valid kmers (kmers which arent in plasmidStr), 
    and tests that they are correct (not in the original)

    Args:
        plasmidStr: the plasmid to read in
        kmer: the kmer size we are looking for
    """
    # Get the possible kmers associated with the plasmid file
    kmers = KmerUtil.GetKmersNotAppearing(plasmidStr,kmer)
    KmerUtil.AssertValidKmers(plasmidStr,kmers)
    # POST: everything is hunky dory.
    return kmers

def GetEmbossPrimers(mCycle,outDir,maxPsuedoLen=None,**kwargs):
    """
    Given a primer cycle (ordered list of kmers such that overlap between 
    i and i+1 isnt in our primer), gets the primers associated with an emboss
    run. These are (nominally) within the melting temperature we want
    Args:
        Cycle: list of primers 
        maxPsuedoLen: maximum pesueoPrimer Length. if non, we done have one
        outDir : where to put the eprimer32 .fasta output
        **kwargs: fed directly into emboss; see the Scripts.Emboss.EmbossAdapter
    Returns:
        list of all the sequences matching the melting temperatures
        set by the emboss utility 
    """
    psuedoPrimer = "".join(mCycle)
    if (maxPsuedoLen is not None and
        len(psuedoPrimer) > maxPsuedoLen):
        psuedoPrimer = psuedoPrimer[:maxPsuedoLen]
    primerFile = outDir + "Primer.fasta"
    embossFile = outDir+"ExamplePrimers.fasta"
    allPrimers=EmbossUtil.GetPrimersFromEmboss(psuedoPrimer,primerFile,
                                               embossFile,**kwargs)
    # Call emboss, get the primers of the size/melting temp we want
    # XXX should make these parameters.
    fwds = [p.forward for p in allPrimers]
    back = [p.backwards for p in allPrimers]
    allPrimers = fwds + back
    # get just the unique ones, by sequence
    unique = list(set(allPrimers))
    return unique

def ReadFileAsPlasmid(inputFile):
    """
    Treats the input string as a file to read in, returns it
    """
    return EmbossUtil.ReadSimpleSeqFile(inputFile)

def GetPlasmidAndKmersNotOccuring(inputFile,kmerLen):
    """
    Given an input (sequence) file and a kmerlength, gets the plasmid and all 
    the kmers not appearing in the input file

    Args: 
        inputFile: the file to read in
        kmerLen: the kmers we are looking for
    Returns: 
        tuple of <plasmid,kmers>
    """
    mPlasmid = ReadFileAsPlasmid(inputFile)
    kmers = GetValidKmers(mPlasmid,kmerLen)
    return mPlasmid,kmers

def GetPossibleKmersFromGraph(mGraph,desiredPrimerLen):
    """
    Given a graph from GraphUtil.GetKmerGraph, gets possible 
    primers of length 'desiredPrimerLen'. Doesn't get them all; finds all
    shortest paths through the graph

    Args:
        mGraph: each node is a kmer, each edge represents the overlap between
        two nodes is also a kmer (ie: also valid / orthogonal to plasmid /
        ok for overhand). Output of GraphUtil.GetKmerGraph
        desiredPrimerLen: length of the desired primers
    Returns:
        a list of the possible primers
    """
    # get the size of the kmer from the first edge, first node
    edges = mGraph.edges()
    kmerLen = len(edges[0][0])
    # if we need primers of length <= 2*node side, just return all edges
    if (desiredPrimerLen <= 2*kmerLen):
        allPrimers = [ KmerUtil.Kmers(e[0]+e[1],desiredPrimerLen)
                       for e in edges]
        # return the flattened list
        Primers = [ Primer(k,None) for listV in allPrimers for k in listV]
    else:
        # POST: more work to do
        # use the GraphUtil method to enumerate all simple paths between
        # all pairs of nodes,then get what we need from those
        Primers = [Primer(k,None) for k in
                   GraphUtil.GetValidPrimersFromGraph(mGraph,desiredPrimerLen)]
    return Primers

def GetOverhangingPrimersGraphAndCycle(inputFile,primerLen,kmerLen):
    """
    Given a plasmid, gets everything we need to call emboss
    Args:
        inputFile: where the simple primer comes from
        primerLen: length of the primers we want
        kmerLen: length of the kmers we want
    returns:
        tuple of <plasmidStr,kmersValid,kmerGraph,PossiblePrimers>
    """
    mPlasmid,kmers = GetPlasmidAndKmersNotOccuring(inputFile,kmerLen)
    # get the graph corresponding to the kmer path overlap
    mGraph = GraphUtil.GetKmerGraph(kmers)
    primers = GetPossibleKmersFromGraph(mGraph,primerLen)
    return mPlasmid,kmers,mGraph,primers
        
def GenerateOverhangingPrimers(inputFile,primerLen,kmerLen):
    """
    Given an input plasmid sequence and primer specifications, generates the
    overhanging primers we could use

    Whatever it returns, makes sure that none of the kmers overlap, or throws
    an erorr.

    Args:
        inputFile: the input file path, assumed just a simple plasmid file (ie:
        with a raw string for input)

        kmerLen: the length of the minimum kmer to use. 

        primerLen: the length of the primer we want

    Returns:
        A 'OverhangPrimerInfo' Object
    """
    mPlasmid,kmers,mGraph,mPrimers = \
        GetOverhangingPrimersGraphAndCycle(inputFile,primerLen,kmerLen)
    # check all our primers are OK.
    KmerUtil.AssertPrimersValid(mPlasmid,kmerLen,mPrimers)
    # POST: all OK
    toRet = OverhangPrimerInfo(mPlasmid,kmers,mGraph,mPrimers)
    return toRet

def GetSequencesWithMinimumDimerization(possibleSeqs,fudge=0):
    """
    Gets the sequences in "possibleSeqs" with the minimum dimerization, accd
    to IDTs method (no penalty for gap extenson)
    
    Args:
        possibleSeqs : list of possible sequences
        fudge : pos number. by default, we get scores with the minimum 
        dimerization. if fudge is greater than 0, gets all scores within 'fudge
        of the min
    Returns:
        list of scores and sequences meeting the recquirements
    """
    sortedSeq,sortedScores = GetSortedScores(possibleSeqs)
    # only get the minimum score
    minV = min(sortedScores)
    # get where the sequences have the minimum homo dimerization
    minIdx = [i for i in range(len(sortedSeq))
              if sortedScores[i] <= minV + fudge]
    minSeq = [sortedSeq[i] for i in minIdx]
    minScore = [sortedScores[i] for i in minIdx]
    return minSeq,minScore

def GetSequenceAndSortByMeltingTemp(Sequences,DesiredMelt):
    """
    Loops through Sequences (primer objects) and sorted by difference from 
    the desired melting temperature, using IDT's method of calculating primers

    Args:
        Sequences: list of primers
        deisredMelt: the desired melting temperature, celsius
    Returns:
        List of Sequences, sorted by distance**2 from DesiredMelt
    """
    # now get all the scores
    meltingTemps = []
    for i,s in enumerate(Sequences):
        temperature = MeltingUtil.GetIdtMeltingTemperatureForPCR(str(s))
        s.SetTemperature(temperature)
        meltingTemps.append(s.temp)
    # sort by the mean...
    sortIdx = np.argsort( (np.array(meltingTemps)-DesiredMelt)**2 )
    sequenceSorted = [Sequences[i] for i in sortIdx]
    return sequenceSorted

def GetMink(mPlasmid,thresh=50,CheckCircular=True):
    """
    Gets the minimum k with at least thresh kmers not occuring in Plasmid
    
    Args:
        mPlasmid: plasmid as a  string
        thresh: the minimum number of kmers to find, for whatever k we return
        CheckCircular: if true, checks the kmers in the circular and overhang
        plasmid
    """
    return KmerUtil.SmallestNonCoveredKmer(mPlasmid,thresh=thresh,
                                           CheckCircular=CheckCircular)


def GetBestOverhangs(inputFile,primerLen,DesiredMelt,kmerLen=None,fudge=0):
    """
    Gets the 'best' overhangs, where best is lowest self-dimerization, no
    'kmerLen'-mers in common with the plasmid, and close to DesiredMelt.

    Uses GenerateOverhangingPrimers, so makes sure none of the kmers it uses
    overlap, or throws and error

    Args:
        inputFile: see GenerateOverhangingPrimers
        primerLen: see GenerateOverhangingPrimers
        desiredMelt: desired melting temperature
        kmerLen: see GenerateOverhangingPrimers. If none, figures out
        the minimum kmer which is reasonable, by 'GetMinkK'
        fudge: see GetSequencesWithMinimumDimerization
    Returns:
        OverhangPrimerInfo object, '.AllPrimers' is a sorted list of 
        all possible, '.BestPrimers' are those with low self-dimerization
    """
    if (kmerLen == None):
        kmerLen = GetMink(ReadFileAsPlasmid(inputFile))
    # get the  information object (this checks that all the primers we return
    # are OK )
    info = GenerateOverhangingPrimers(inputFile,primerLen,kmerLen)
    # now we filter by dimerization score...
    minSeq,scores = GetSequencesWithMinimumDimerization(info.AllPrimers,
                                                        fudge=fudge)
    for seq,score in zip(minSeq,scores):
        seq.SetAlignment(score)
    # and by melting temperature
    sortedSeq = GetSequenceAndSortByMeltingTemp(minSeq,DesiredMelt)
    info.SetBestPrimers(sortedSeq)
    return info
    
def GetSortedScores(sequences):
    """
    Given a list of sequences (e.g. primers), gets the sorted
    list of their self-dimerization scores, low to high
    Args:
        sequences: the list of sequence
    """
    scores = [AlignUtil.GetBestSelfDimerAlignmentScore(str(s))
              for s in sequences]
    bestDimerScores = np.array(scores)
    # sort the dimer score to get the best possible (worst dimer alignment)
    sortIdx = [i for i in np.argsort(bestDimerScores) if scores[i] is not None]
    # print them all out
    sortedSeq = [sequences[i] for i in sortIdx]
    sortedScores = bestDimerScores[sortIdx]
    return sortedSeq,sortedScores

def SafeKmerToAdd(info,kToAdd=None,MaxRepeats=2,MaxBadKmers=0):
    """
    Given an info object (from GetBestOverhangs), gets the best kToAdd-mer to 
    add to the (5' end) of the best sequence (first seq in info.BestPrimers). 
    This assumes we are only adding to get to IDTs minimum, doesnt look
    at what happens to the melting temperature

    Throws an error if it cant find a safe kmer.

    Args:
        info: the information object from GetBestOverhangs
        kToAdd : the kmer to add. If none, then defaults to 15-actualK

        MaxRepeats: maximum number of allowable repeats (defaults to 2, so 
        AA is ok in whatever we add, but not AAA)

        MaxBadKmers: maximum bad kmers, non negative number. If 0, no 
        kmers in the resulting 'concatenated' kmer will be in the plasmid.
        If 1, at most one overlapping kmer
    Returns:
        tuple of <kmer to add,number of times kmer occurs in plasmid>
    """
    bestPrimer = info.BestPrimers[0]
    # get the lowest 3-mer (so we can send to IDT)
    kmerAlpha = info.PossibleKmers
    if (kToAdd is None):
        # determine from IDTs default of 15
        fullPrimerLen = IDT_MIN_NUM
        # determine kToAdd based on our alphabet...
        kToAdd = fullPrimerLen - len(bestPrimer)
    bestKmersSorted = KmerUtil.LeastOccuringKmers(info.Plasmid,kToAdd)
    # need to know how big valid kmers are
    kmerSize = info.AlphabetK
    kmerStr = str(bestPrimer)
    badKmers = set(["".join([a for _ in range(MaxRepeats+1)])
                    for a in KmerUtil.DEF_BASES_DNA])
    # assume we dont find it by default
    Found=False
    for k,score in bestKmersSorted:
        if k in badKmers:
            continue
        concatPrimer = k + kmerStr
        kmersTmp = set(KmerUtil.Kmers(concatPrimer,kmerSize))
        if len(kmersTmp - kmerAlpha) <= MaxBadKmers:
            # then adding this kmer is OK, all of its kmers are OK
            Found = True
            break
    assert Found , "Couldn't find kmer to add to make..."
    # make sure we didnt add too much crap...
    overlap = (set(KmerUtil.Kmers(k+kmerStr,kmerSize)) - kmerAlpha)
    assert len(overlap) <= MaxBadKmers
    return k,score


def GetPrimersSortedByAlignmentScore(mInfo):
    """
    Given output of GenerateOverhangingPrimers, gets the sorted primers 
    and scores associated with self dimerization
    Args:
        mInfo: primer info object
    Returns:
        tuple of <primerList,scoreList>, where index i refers to the i
        *lowest* scoring primer
    """
    # get the worst-possible dimerization score (best dimer alignment)
    # for each sequence
    allSeq = mInfo.AllPrimers
    return GetSortedScores(allSeq)


def ConcatToPrimers(PrimerStr,ForwardPrimer,ReversePrimer,
                    addSpacer=False,addDbcoAndBio=False,
                    **kwargs):
    """
    see ConcatPrimerTo1607Fand3520R for arguments

    XXX todo: docs
    """
    revComp = KmerUtil.ReverseComplement(PrimerStr)
    if (addSpacer or addDbcoAndBio):
        # convert the strings to lists
        func = lambda x: list(x)
    else:
        func = lambda x : x
    if (not addSpacer):
        # nothing funny, just concat the strings
        add = ""
    else:
        add = [IdtUtil.IdSpacer()]
    # do we need to add the dbco and biotin?
    if (not addDbcoAndBio):
        AddFwdPre = ""
        AddRevPre = ""
    else:
        AddFwdPre = [IdtUtil.Dbco5Prime()]
        AddRevPre = [IdtUtil.Biotin5Prime()]
    fwd = func(AddFwdPre) + func(PrimerStr) + func(add) + func(ForwardPrimer)
    rev = func(AddRevPre) + func(revComp) + func(add) + func(ReversePrimer)
    # may need to ignore spacers, etc....
    offset = len(AddFwdPre)
    # check and make sure we didn't screw anything up...
    AssertPrimersAndOverhangCorrect(fwd[offset:],rev[offset:],
                                    PrimerStr,ForwardPrimer,ReversePrimer)
    return fwd,rev
    

def ConcatPrimerTo1607Fand3520R(primerStr,
                                addSpacer=False,
                                addDbcoAndBio=False,**kwargs):
    """
    Concatenates the given primer to the 5' end of the 1607F and the 5' end 
    of the 3520R 

    Asserts that the overhang and its reverse compement are OK, as well as the
    actual primers werent screwed up.
    
    Args:
       primerStr: the 5' to 3' primer to add
       addSpacer: if True, returns string as a list, with an idSpacer separating
       AddDbcoAndBiotin: if True, returns string as a list, with 5' DBCO
       on the forward, and 5' Biotin on the reverse

       **kwargs: bassed to CommonPrimers.Get1607Fand3520R
    returns:
       The forward and reverse primer
    """
    pFwd,pRev = CommonPrimerUtil.Get1607FAnd3520R(**kwargs)
    return ConcatToPrimers(primerStr,pFwd,pRev,addSpacer,addDbcoAndBio,
                           **kwargs)

def AssertPrimersAndOverhangCorrect(fwdToCheck,revToCheck,ovh,fwd,rev):
    """
    Given an actual overhang (small), fwd and reverse primer (to concatenate
    with the overhang, checks that

    (1) the overhang and its reverse complement are in their right place
    in the forward and reverse primer (we assume overhang goes on 5')
    (2) same with the fwd and reverse (we assume overhang goes on 5')

    Args:
        fwdToCheck: the forward primer we want to make; this is what we check
        revToCheck: ibid, bu reverse
        ovh: the overhang we want (using this as 'gold standard')
        fwd: the forward primer we want (using this as 'gold standard')
        rev: the reverse primer we want (using this as 'gold standard')
    """
    nFwd = len(fwd)
    nRev = len(rev)
    # in case we have complicated stuff (e.g. biotin, whatever), go ahead
    # and join the strings together
    sanit = lambda x: "".join(x)
    # check the primers themselves    
    assert sanit(fwdToCheck[-nFwd:]) == sanit(fwd)
    assert sanit(revToCheck[-nRev:]) == sanit(rev)
    # check the overhanging revions for the primer and its reverse complement
    lenOvh = len(ovh)
    revComp = KmerUtil.ReverseComplement(ovh)
    assert sanit(fwdToCheck[:lenOvh]) == sanit(ovh)
    assert sanit(revToCheck[:lenOvh]) == sanit(revComp)

def ConcatAndSave(mPrimer,baseDir,Name,
                  ForwardSequence,ReverseSequence,
                  addSpacer=False,addDbcoAndBio=False,
                  **kwargs):
    fwdSpacer,revSpacer = ConcatToPrimers(mPrimer,ForwardSequence,
                                          ReverseSequence,
                                          addSpacer,addDbcoAndBio,**kwargs)
    # make orders, including the spaces
    spacerStr = "_Spacer" if addSpacer else ""
    spacerStr += ("_DBCO_BIO" if addDbcoAndBio else "")
    mOrders = [[fwdSpacer,"F_" + Name + spacerStr],
               [revSpacer,"R_" + Name + spacerStr]]
    # must have 100nm scale for IdSp as of 4/1/2016
    scales = IdtUtil.Scales
    scale = scales._100NM if (addSpacer) else scales._25NM
    spacerOrders = IdtUtil.SequencesAndNamesTuplesToOrder(mOrders,
                                                          Scale=scale)
    IdtUtil.PrintAndSave(spacerOrders,"./" + Name + spacerStr + ".txt")


def ConcatTo1607FAnd3520RAndSave(mPrimer,baseDir,Name,addSpacer=False,
                                 addDbcoAndBio=False,**kwargs):
    """
    Given a primer, concatensates it to the 1607F and 3520R primer, optionally
    adds spacers and labels, and 

    Args:
        mPrimer: primer, as simple string, to use
        baseDir: where 1607F and 3520R live
        Name: base name for the primer; will be tacked into 1607F and 3520R
        as appropriate, with additions for spacers and dbco.biotin
    
        addSpacer:  if true, adds a spacer between the primer and the 
        'normal' 1607F and 3520R
        addDbcoAndBio: if True, adds a 5' dbco and 5' biotin to the forward
        and reverse primers,respectively
    Returns: nothing
    """
    pFwd,pRev = CommonPrimerUtil.Get1607FAnd3520R(base=baseDir,
                                                  **kwargs)
    ConcatAndSave(mPrimer,baseDir,Name,
                  pFwd,pRev,addSpacer,addDbcoAndBio,
                  **kwargs)
    
def CreateOverhangsFor1607F(inputFile,baseDir,desiredPrimerLen,desiredMeltTemp,
                            Name,MakeSpacerFile=False,MakeLabelledFile=False,
                            CheckAfterChoosing=True,MakePrimerFile=True,
                            ChooseFunc = lambda x: x[0],
                            MaxPlasmidAlignment=None,**kwargs):
    """
    Args:
        inputFile: where the primer file lives
        BaseDir: see ConcatTo1607FAnd3520RAndSave
        Name: see ConcatTo1607FAnd3520RAndSave
        desiredMeltTemp : what melting temperature we want; tries to 
        get as close as the other constraints allow

        MakeSpacerFile: if true, makes a file with spacers bewteen, see 
        ConcatTo1607FAnd3520RAndSave

        MakeLabelledFile: if true, makes a file with dbco/bio labels, see 
        ConcatTo1607FAnd3520RAndSave

        ChooseFunc: Given a list of primers sorted by absolute dist from melting
        temperature, this returns a single primer to use/save, based on 
        other criteria (e.g. GC content, etc). Default to just closest by melt

        CheckAfterChoosing: If true, checks that the primer returned is still
        OK (ie: no kmer-overlaps with the plasmid) after we return

        MaxPlasmidAlignment: maximum idt alignment score between a valid
        overhanging primer and the plasmid. If none, defaults to the kmer size
    
        kwargs: passed to GetBestOverhangs
    """
    # determine the best overhang (low self-dimerization, close to desired
    # melting temperature)
    inf = GetBestOverhangs(inputFile,desiredPrimerLen,desiredMeltTemp,
                           **kwargs)
    # get the actual best primers
    PrimersRaw = inf.BestPrimers
    if (MaxPlasmidAlignment is None):
        MaxPlasmidAlignment = inf.AlphabetK
    # XXX TODO check global alignments with entire plasmid.
    """
    InStr = EmbossUtil.ReadSimpleSeqFile(inputFile)
    AlignFunc = AlignUtil.GetIdtAlignments
    print("acqui...")
    print("n={:d}".format(len(PrimersRaw)))
    Scores = [AlignFunc(InStr,str(p),one_alignment_only=True)
              for p in PrimersRaw]
    print("alla...")
    primers = [p for i,p in enumerate(PrimerRaw)
               if Scores[i] <= MaxPlasmidAlignment]
    """
    primers = PrimersRaw
    # POST: all overhangs align poorly with the plasmid. Yay!
    bestPrimer = ChooseFunc(primers)
    # re-calculate the alignment and temperature, in case we did something
    # funky
    bestPrimer.SetAlignment()
    PrimerTemp = MeltingUtil.GetIdtMeltingTemperatureForPCR(str(bestPrimer))
    bestPrimer.SetTemperature(PrimerTemp)
    print("We Picked (Primer, Melting temperature, Alignment)")
    print(bestPrimer.seq,bestPrimer.temp,bestPrimer.SelfAlignment)
    if(CheckAfterChoosing):
        KmerUtil.AssertPrimersValid(inf.Plasmid,inf.AlphabetK,[bestPrimer])
    idtStr = str(bestPrimer)
    if MakePrimerFile:
        ConcatTo1607FAnd3520RAndSave(idtStr,baseDir,addSpacer=False,Name=Name)
    if (MakeSpacerFile):
        ConcatTo1607FAnd3520RAndSave(idtStr,baseDir,addSpacer=True,Name=Name)
    if (MakeLabelledFile):
        ConcatTo1607FAnd3520RAndSave(idtStr,baseDir,addSpacer=True,
                                     addDbcoAndBio=True,
                                     Name=Name)
    inf.SetAllPossiblePrimers(primers)
    return inf,bestPrimer

def ChooseFirstWithout3PrimeGC(primers):
    mList = [p for p in primers if (p[-1].upper()) not in ["G","C"]]
    return mList[0]

def ChooseFirstPrimerMatching(primers,ToMatch,
                              AddTo3Prime="",
                              AddTo5Prime=""):
    mList = [p for p in primers if (ToMatch.lower() in str(p).lower())]
    toRet = mList[0]
    toRet.seq = AddTo5Prime + toRet.seq + AddTo3Prime
    return toRet

def ChooseOverhangByBothEnds(primers):
    """
    Chooses overhang to have a GC-rich 5' end, and a GC-poor 3' end 
    (endsureing 

    Args:
        primers: list of primers.
    Returns:
        first primer in list matching what we want 
    """
    Size = 5
    MaxGCThreePrime = 1
    MinGCFivePrime = 3
    # Get the GC Content
    CountFunc = lambda x: x.upper().count("C") + x.upper().count("G")
    ThreePrimeEnds = [CountFunc(p[-Size:]) for p in primers]
    FivePrimeEnds = [CountFunc(p[:Size]) for p in primers]
    # get only those meeting the recquirements
    mList = [p for i,p in enumerate(primers) if
             (ThreePrimeEnds[i] <= MaxGCThreePrime)
             and
             (FivePrimeEnds[i] >= MinGCFivePrime)]
    return mList[0]

def GetBestAlignmentsWithPlasmid(inf):
    """
    Given return of CreateOverhangsFor1607F, gets *all* the alignment scores
    associated with the plasmid
    """
    scores = []
    Plasmid = inf.Plasmid
    for i,primer in enumerate(inf.AllPrimers):
        score = AlignUtil.GetEbiAlignments(str(primer),Plasmid,
                                           one_alignment_only=True)[0].score
        scores.append(score)
    return scores


def GetPrimersWithoutOverhangs(Plasmid,ProductLength,
                               SliceOther,OtherIsReverse,PrimerLength,Overhang,
                               Name):
    """
    Args: see CreatePrimer
    """
    # overhang is length + *both* abasic site (the overhangs will anneal like
    # 5'                      -- overhang -- abasic -- Forward - 3'
    # 3' -- Reverse -- abasic -- overhang -- 5'
    # so the abasic sites add an effective TWO base pairs
    OverhangLength = len(Overhang) + 2
    # need zero based to get thi
    StartOther = SliceOther.start
    EndOther = SliceOther.stop
    LengthOtherPrimer = abs(EndOther-StartOther)+1
    TotalPrimerLength = PrimerLength +LengthOtherPrimer +OverhangLength
    # get the remaining length
    DistanceBetweenPrimers = ProductLength-TotalPrimerLength
    # figure out where *this* primer should start
    if (OtherIsReverse):
        # get the start
        StartOfOther = min(StartOther,EndOther)
        Start = StartOfOther - DistanceBetweenPrimers - PrimerLength
    else:
        Start = StartOther+ LengthOtherPrimer + DistanceBetweenPrimers
    End = Start + PrimerLength
    PrimerSlice = slice(Start,End,1)
    # determine which slice goes where
    if (OtherIsReverse):
        ReverseSlice = SliceOther
        ForwardSlice = PrimerSlice
    else:
        ReverseSlice = PrimerSlice
        ForwardSlice = SliceOther
    # cool, now just print out the sequences
    PrimerForwardSeq = Plasmid[ForwardSlice]
    PrimerReverseSeq = KmerUtil.ReverseComplement(Plasmid[ReverseSlice])
    return PrimerForwardSeq,PrimerReverseSeq,ForwardSlice,ReverseSlice

def CreatePrimer(Plasmid,ProductLength,
                 SliceOther,OtherIsReverse,PrimerLength,Overhang,Name):
    """
    Given an overhang and desired length, creates the primers needed
    (full 'suite': normal, overhang, etc

    Args:
        Plamid: full plasmid we want
        Product Length: how long the product should be, *including* the 
        abasic site and overhang
        SliceOther: indices for the other primer
        OtherIsReverse: if true, other is reverse complement of its slice
        Otherwise, this primer is reverse complement
        PrimerLength: how long the (new) primer should be, without the 
        overhang
        Overhang: sequence of the overhang
        Name: base name for this construct
    """
    PrimerForwardSeq,PrimerReverseSeq,ForwardSlice,ReverseSlice = \
        GetPrimersWithoutOverhangs(Plasmid,ProductLength,
                                   SliceOther,OtherIsReverse,PrimerLength,
                                   Overhang,Name)
    # for without the overhang, dont add *anything*
    ConcatAndSave("",baseDir="./",Name=Name,
                  ForwardSequence=PrimerForwardSeq,
                  ReverseSequence=PrimerReverseSeq,
                  addSpacer=False,addDbcoAndBio=False)
    # add the overhangs and labels
    ConcatAndSave(Overhang,baseDir="./",Name=Name,
                  ForwardSequence=PrimerForwardSeq,
                  ReverseSequence=PrimerReverseSeq,
                  addSpacer=True,addDbcoAndBio=False)
    ConcatAndSave(Overhang,baseDir="./",Name=Name,
                  ForwardSequence=PrimerForwardSeq,
                  ReverseSequence=PrimerReverseSeq,
                  addSpacer=True,addDbcoAndBio=True)


             
        
