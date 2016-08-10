import Bio.SubsMat.MatrixInfo
from Bio import pairwise2
import KmerUtil
from Bio.pairwise2 import format_alignment

DEF_MATCH = 1
DEF_MISMATCH = -1
DEF_GAP = -1
 # idt evidentally uses -2 for opening a gap or mismatching
IDT_DEF_GAP_MIS = -2
# Parameters from European Bioinformatics Institute,
# http://www.ebi.ac.uk/Tools/psa/emboss_needle/nucleotide.html
# 2016-8-9:
EBI_GAP_OPEN = -10
EBI_GAP_EXTEND = -0.5
# see : http://osdir.com/ml/science.biology.emboss/2005-12/msg00028.html
EBI_MISMATCH = -4
EBI_MATCH = 5

class AlignmentInfo:
    def __init__(self,s1,s2,score,startIdx,endIdx):
        """
        Wrapper class for Biopython. just wraps the sequences we need.
        All args come directly from pairwise2.align.globalXX
        
        Unit tested by TestUtil.TestAlignments.TestReverseComplementAlignments
        Args:
            s1: alignment of s1
            s2: alignment of s2
            scoore: score of the alignment
            startIndex: start index, 4th arg returned by pairwise
            endIdx: start index, 5th arg returned by pairwise
        """
        self.s1 = s1
        self.s2 = s2
        self.score = score
        self.startIdx = startIdx
        self.endIdx = endIdx
    def AlignmentTuple(self):
        """
        Return the tuple which pairwise2.align.global returns. Useful for
        (e.g.) pretty printing
        """
        return (self.s1,self.s2,self.score,self.startIdx,self.endIdx)
    def __str__(self):
        return format_alignment(*self.AlignmentTuple())
    def __repr__(self):
        return str(self)

def Sanitize(Seq):
    """
    Sanitizes (strips out trailing/starting whitespace, lowercase) 
    a given sequence

    Args:
        Seq: The sequence to sanitize
    Returns:
        The sanitized sequence
    """
    return Seq.strip().upper()

def GetIdtAlignments(Seq1,Seq2,MismatchScore=IDT_DEF_GAP_MIS,
                     GapOpen=IDT_DEF_GAP_MIS,
                     GapExtend=0,**kwargs):
    """
    Gets the alignment scores for two sequences,using (by default) IDT's params,
    
    Args:
        Seq1,Seq2: align Seq1 to Seq2. *both should be DNA
        Other args: see GetBestSelfDimerAlignmentScore
    Returns:
        maximum over all alignment scores
    """
    alignments = AlignmentScores(Seq1,Seq2,MismatchScore=MismatchScore,
                                 GapOpen=GapOpen,GapExtend=GapExtend,**kwargs)
    return alignments


def GetEbiAlignments(Seq1,Seq2):
    """
    Gets the EBI (European Bioinformatics Institute) local alignment on DNA,
    using defaults listed on 
    ebi.ac.uk/Tools/psa/emboss_needle/help/index-nucleotide.html
    
    Args: 
        See AlignmentScores: both are DNA
    """
    return AlignmentScores(Seq1,Seq2,
                           MatchScore=EBI_MATCH,
                           MismatchScore=EBI_MISMATCH,
                           GapOpen=EBI_GAP_OPEN,
                           GapExtend=EBI_GAP_EXTEND)


def GetBestSelfDimerAlignmentScore(Seq,MismatchScore=IDT_DEF_GAP_MIS,
                                   GapOpen=IDT_DEF_GAP_MIS,
                                   GapExtend=0,**kwargs):
    """
    Gets the best (highest) self-dimer alignment for the given sequence with
    its reverse complement (this states 'how likely is the sequence to 
    bind to itself). 

    By default, similar to what the Homo-Dimer Analysis of Idt 
    (http://www.idtdna.com/calc/analyzer , look for "Self-Dimer") does. It 
    peanlizes any gas by two 

    e.g. TAGGACCACTCG -> 2 are most according to ids

    Unit tested by TestUtil.TestAlignments.TestReverseComplementAlignments
    
    Args:
        Seq: Sequence to align with its reverse
        Others: see AlignmentScores. Note default has no penalties
    Returns:
        score from alignment. If using default arguments, this is number 
        of base-pair matches, less 2 for the start of any gap
    """
    alignments = AlignSelfWithReverseComplement(Seq,
                                                MismatchScore=MismatchScore,
                                                GapOpen=GapOpen,
                                                GapExtend=GapExtend,
                                                **kwargs)
    bestScore = max([a.score for a in alignments])
    return bestScore
        
def AlignSelfWithReverseComplement(Seq,MismatchScore=IDT_DEF_GAP_MIS,
                                   GapOpen=IDT_DEF_GAP_MIS,GapExtend=0,
                                   **kwargs):
    """
    Gets an alignment score for the sequence with itself reversed, complemented.

    Unit tested implicitly by
    TestUtil.TestAlignments.TestReverseComplementAlignment

    Args:
        Seq: Sequence to align with itself
        Others: See AlignmentScores
    Returns:
        List of possible alignments as AlignmentInfo objects
    """
    ReverseComp = KmerUtil.ReverseComplement(Seq)
    return AlignmentScores(Seq,ReverseComp,MismatchScore=MismatchScore,
                           GapOpen=GapOpen,GapExtend=GapExtend,**kwargs)
                           

def AlignmentScores(Seq1,Seq2,MatchScore=DEF_MATCH,MismatchScore=DEF_MISMATCH,
                    GapOpen=DEF_GAP,GapExtend=DEF_GAP,SanitizeSeqs=True,
                    **kwargs):
    """
    Align two sequences locally. 

    Unit tested implicitly by
    TestUtil.TestAlignments.TestReverseComplementAlignment

    Args:
        Seq1: First Sequence
        Seq2: Second Sequence
        MatchScore: Amount to add per match
        MismatchScore: Amount to add per mismatch. Usually <0
        GapOpen : Amount to add per gap Open. Usually <0
        GapExtend: Amount to add per gap extension (given open). Usually <0
        SanitizeSeqs: If true, calls the sanitize function on input strings
        **kwargs: passed to localms
    Returns:
        List of possible Alignments as AlignmentInfo objects
    """
    # see http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
    # look for 'globalms'
    if (SanitizeSeqs):
        Seq1 = Sanitize(Seq1)
        Seq2 = Sanitize(Seq2)
    alignments =  pairwise2.align.localms(Seq1,Seq2,MatchScore,MismatchScore,
                                          GapOpen,GapExtend,**kwargs)
    if (len(alignments) == 0):
        # no alignment possible
        return [AlignmentInfo(Seq1,Seq2,None,None,None)]
    else:
        return [AlignmentInfo(*a) for a in alignments]

                                             
