import KmerUtil

class Scales:
    _250NM = "250nm"
    _100NM = "100nm"
    _25NM = "25nm"

class Purifications:
    HPLC = "HPLC"
    NONE = "STD"
    PAGE = "PAGE"

class IdtOrder:
    def __init__(self,ProductName,Sequence,Scale=Scales._25NM,
                 Purification=Purifications.NONE):
        """
        Records all the information needed to specify an IDT order
        Args:
            ProductName: name of the product
            Sequence: the actual sequence, typically a string (if nothing fancy,
            such as "ATCG") or list of characters (if fancy, like with biotin:
            ["/Biotin/","A",T"])
            
            yieldScale: valid yield to use. From class Yields
            purification: purification to use. From class Purifications. 
            If DBCO or Biotin is present, overrides
        """
        self.name = ProductName
        self.seq = Sequence
        self.scale = Scale
        self.purification = Purification
        # biotin or DBCO recquires HPLC (not even PAGE works)
        if ( (Biotin5Prime() in self.seq) or (Dbco5Prime() in self.seq)):
            self.purification = Purifications.HPLC
    def __str__(self):
        return"{:s}\t{:s}\t{:s}\t{:s}".format(self.name,
                                              IdtSeqStr(self.seq,
                                                        FivePrimeMarkers=False),
                                              self.scale,
                                              self.purification)
    def __getitem__(self,index):
        return self.seq[index]
    def __len__(self):
        return len(self.seq)
    
def PrintSequencesToFile(seqs,outfile,CommonAssertions=True):
    """
    Prints a series of sequences to files

    Args:
       seqs: list of sequences (class IdtOrder) to print
       outfile: where to print the sequences
       CommonAssertions: if true, throws error in (probably) unwanted situations
       (ie: two sequences with the same content or name)
    Returns: 
       nothing
    """
    if (CommonAssertions):
        # check that the names and sequences are unique
        seqNames = [s.name for s in seqs]
        seqContent = ["".join(s.seq) for s in seqs]
        # dont want any duplicate names or sequences
        assert len(seqNames) == len(set(seqNames)) , "Identical Sequence Names"
        assert len(seqContent) == len(set(seqContent)) , "Identical Sequences"
    with open(outfile,'w') as f:
        mStrs = [str(s) for s in seqs]
        f.write("\n".join(mStrs))

def IdtSeqStr(Seq,FivePrimeMarkers=True):
    """
    Converts a sequence (e.g. simple string or list of elements)
    into IDT pretty printing format
    
    Args:
        Seq: The sequene to print. Each element is considered a letter.
        If using labels, just make a list, each list element is a string
        (e.g. ["A","/Biotin/","T"] or "ATC" are acceptable, but not 
        "A/Biotin/T")
    """
    # group by threes
    seqLen = len(Seq)
    byThrees = [ "".join(Seq[i:min(seqLen,i+3)]) for i in range(0,seqLen,3)]
    # print off with spaces, and 5/ 3/ markers
    joinedIdt =" ".join(str(t) for t in byThrees)
    if (FivePrimeMarkers):
        mStr = "/5/" + joinedIdt + "/3/"
    else:
        mStr = joinedIdt
    return mStr

def SequencesAndNamesToOrder(Seqs,Names,**kwargs):
    """
    Converts the given names and sequenes to an order

    Args:
        Seqs: sequences, as string
        Names: Names, as strings
        **kwargs: see SequenceAndReverseComplementOrders
    Returns: list of orders
    """
    assert (len(Seqs) == len(Names))
    toRet = []
    for seq,name in zip(Seqs,Names):
        toRet.append(IdtOrder(name,seq,**kwargs))
    return toRet

def SequencesAndNamesTuplesToOrder(SeqsNames,**kwargs):
    """
    Converts the given names and sequenes to an order

    Args:
        SeqsNames: list of tuples, each is (sequence, name)
        **kwargs: see SequenceAndReverseComplementOrders, passed to order
    Returns: 
        slist of orders
    """
    seqs = [s[0] for s in SeqsNames]
    names = [s[1] for s in SeqsNames]
    return SequencesAndNamesToOrder(seqs,names,**kwargs)



def SequenceAndReverseComplementOrders(Seq,NameFwd,NameRev,**kwargs):
    """
    Gets the sequence and its reverse complement as orders

    Args:
        Seq: Sequence to use
        NameFwd: Name of the forward order
        NameFwd: Name of the reverse order
        **kwargs: passed directly to the order, assumed the same
    Returns:
        list of orders
    """
    revComp = KmerUtil.ReverseComplement(Seq)
    mSeqs =  [Seq,revComp]
    mNames = [NameFwd,NameRev]
    return SequencesAndNamesToOrder(mSeqs,mNames,**kwargs)

def IdSpacer():
    """
    See: idtdna.com/site/Catalog/Modifications/Category/2

    Return: the string for an abasic spacer. 
    """
    return "/IdSp/"

def Biotin5Prime():
    """
    See: idtdna.com/site/Catalog/Modifications/Category/2

    Returns: the string for a 5' Biotin, with TEG spacer
    """
    return "/5BiotinTEG/"

def Biotin3Prime():
    """
    See: idtdna.com/site/Catalog/Modifications/Category/2

    Returns: the string for a 3' Biotin, with TEG spacer
    """
    return "/3BioTEG/"
    

def Dbco5Prime():
    """
    Returns: the string for a 5' Biotin
    """
    return "/5DbcoTEG/"

def Cy5_dye_3_prime():
    """
    Returns: the string for a 3' Cy3
    """
    return "/3Cy5Sp/"

def atto_633():
    """
    Returns: the string for a 3' atto 633 -nhs
    """
    return "/3ATTO633N/"



def PrintOrders(orders):
    """
    Prints the given orders, one on each line

    Args:
        orders: list of IdtOrder objects.
    """
    for o in orders:
        print(o)

def PrintAndSave(orders,fileV):
    """
    Prints and saves the orders we want.

    Args:
        orders:  see PrintOrders:
        fileV: where to save the file
    """
    PrintOrders(orders)
    PrintSequencesToFile(orders,fileV)


def AddDBCOAndBiotin(Sequence):
    """
    Adds a 5' DBCO and a 3' biotin  to the given sequence

    Args:
         Sequence: what to add to
    Returns:
         Sequence as a list of characters, with a DBCO and Biotin appended
    """
    SeqByChar = [l for l in Sequence]
    return [Dbco5Prime()] + SeqByChar + [Biotin3Prime()]
