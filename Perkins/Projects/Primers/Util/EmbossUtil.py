# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
# need to add the utilities class. Want 'home' to be platform independent
from GeneralUtil.python import GenUtilities as pGenUtil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Emboss import Primer3
import Research.Perkins.Projects.Primers.Scripts.Emboss.EmbossAdapter as \
    EmbossAdapter
import AlignUtil

import KmerUtil
import MeltingTemperatureUtil

class Primer:
    def __init__(self,seq,temp=None):
        self.seq= seq
        self.temp = temp
        # alignment score may get set later
        self.SelfAlignment = None
    def SetAlignment(self,a=None):
        """
        Sets the alignment score to a. If 'a' isnt present, uses the
        best self-diemr alignment score method

        Args:
            a: alignment score. If none, auto-calculates
        """
        # sets the alignment score 
        if (a is None):
            a = AlignUtil.GetBestSelfDimerAlignmentScore(self.seq)
        self.SelfAlignment = a
    def SetTemperature(self,t=None):
        if (t is None):
            t = MeltingTemperatureUtil.GetIdtMeltingTemperature(self.seq)
        self.temp = t
    def __eq__(self,other):
        """
        We are equal to other is out sequences are the same
        """
        return other.seq == self.seq
    def __hash__(self):
        """
        Same sequence, same hash
        """
        return hash(self.seq)
    def __len__(self):
        return len(self.seq)
    def __getitem__(self,index):
        return self.seq[index]
    def __str__(self):
        return self.seq
    def __repr__(self):
        return str(self)

class PrimerPair:
    def __init__(self,forwardSeq,forwardTemp,backSeq,backTemp):
        self.forward = Primer(forwardSeq,forwardTemp)
        self.backwards = Primer(backSeq,backTemp)

def GetPrimersFromEmboss(Plasmid,InFile,OutFile,**kwargs):
    """
    Creates a fasta file ('infile') from the input plasmid string, then 
    calls emboss using the **kwargs (see EmbossAdapter)

    Args:
        Primers: simple string of the input plasmid.
        InFile: where to save the (intermediate) fasta file of the plasmid
        Outfile: where to save emboss' output
        **kwargs: see Scripts.EmbossAdapter.RunEmbossAndCreatePrimerFile
    """
    SaveAsFasta(Plasmid,InFile)
    # call the emboss utility with the files
    EmbossAdapter.RunEmbossAndCreatePrimerFile(InFile,OutFile,**kwargs)
    # read back in the primers
    # Read back in the file
    allPrimers = GetPrimersFromFasta(OutFile)
    return allPrimers
    
        
def ReadSimpleSeqFile(fileName):
    """
    Reads a simple sequence file which consists *only* of the sequence.

    Args:
        nVals: list of vertex / connected graph sizes to try
         
        means: mean probability to find a random edge, one per nVals

        stdev: std probability to find a random edge, one per nVals
    """
    with open(fileName) as f:
        plasmidStr = KmerUtil.SanitizeSequence(f.read())
    return plasmidStr

def getRec(id,strV):
    return SeqRecord(Seq(strV.strip(),IUPAC.unambiguous_dna),id=str(id))

def SaveAsFasta(seqs,fileName):
    """
    Saves a list of sequences 'seqs' in the fasta file format
    Args:
       seqs: sequences to save out
       fileName: where to save
    Returns:
       None
    """
    if (type(seqs) is list):
        records =[getRec(id=i,strV=s) for i,s in enumerate(seqs) ]
    else:
        # just a string
        records = getRec(id=0,strV=seqs)
    # POST: records has at least one sequence in it
    with open(pGenUtil.ensureEnds(fileName,'.fasta'), "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

def GetPrimersFromFasta(inputFile):
    """
    Reads in the sequences and melting temperatures of an eprimer32 fasta file

    Args:
       inputFile: where to  look for the file
    Returns:
       list of <PrimerPair> objects, which have the information we need 
    """
    primerList = []
    with open(inputFile) as fileHandle:
        record = Primer3.parse(fileHandle)
        # XXX check is len>0
        primers = record.next().primers
        numPrimers = len(primers)
        size=(numPrimers*2,1)
        seqs= []
        temps = []
        for i,p in enumerate(primers):
            primerList.append(PrimerPair(p.forward_seq,p.reverse_tm,
                                         p.forward_seq,p.reverse_tm))
    return primerList

if __name__ == "__main__":
    generateK = True
    if (generateK):
        generateKmers()
    else:
        inputFile = "./outputDir/myresults.out"
        seqs,temps=getSeqTempFromFasta(inputFile)
        plotEmbossResults(seqs,temps)
