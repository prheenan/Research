# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
import sys

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# #get better names for the common table we will use
# see:
# biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html#Tm_NN
NEAREST_N_ALLAWI_LUCIA_1997= mt.DNA_NN3
# see :
#biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-pysrc.html#salt_correction
SALT_MONO_Owczarzy_2004 = 6
SALT_DIV_Owczarzy_2008 = 7

# tables
DEF_NN_TABLE = NEAREST_N_ALLAWI_LUCIA_1997
DEF_SALT=SALT_MONO_Owczarzy_2004 
# default concentrations
DEF_K =0
# 10uL of 8mM (2mM per dNTP) in 100uL, converted to uM
DEF_dNTPs=10*2/100 * 1000
# all concentrations below in mM, unless otherwise noted
DEF_NA=50
DEF_TRIS=0
# 6uL of 25mM in 100uL
DEF_Mg=(6*25/100)
# oligo concentrations, in nM
DEF_OLIGO = 250
DEF_OTHER = 0
                       
class MeltingParams:
    """
    Class to encapsulate melting parameters
    """
    def __init__(self,Na=DEF_NA,K=DEF_K,Tris=DEF_TRIS,Mg=DEF_Mg,dNTPs=DEF_dNTPs,
                 saltcorr=DEF_SALT,nn_table=DEF_NN_TABLE,dnac1=DEF_OLIGO,
                 dnac2=DEF_OTHER,**kwargs):
        """
        Class to wrap up parameters for melting temperature

        Args:
            Na: Sodioum Concentration [mM]
            K : potassium concentraion [mM]
            Tris: tris concentrtion [mM]
            Mg: magnesium concentation  [mM]
            dNTPs: dNTP concentration [mM]
            saltcorr: salt correction method. See 
            nn_table: the nearest neighbor table
            dnac1: concentration of strand 1 [nM]
            dnac2: concentratin of strand 2 [nM]
            kwargs: passed directly to melting temperature
        """

        self.Na = Na
        self.K = K
        self.Tris = Tris
        self.Mg = Mg
        self.dNTPs = dNTPs
        self.saltcorr = saltcorr
        self.nn_table = nn_table
        self.dnac1 = dnac1
        self.dnac2 = dnac2
        # anything else, also adds
        for key,val in (kwargs.items()):
            setattr(self,key,val)
    def concdict(self):
        """
        Returns all concentrations as a dict
        """
        return dict(Na=self.Na,
                    K=self.K,
                    Tris=self.Tris,
                    Mg=self.Mg,
                    dNTPs=self.dNTPs)

# 3/16/2015:
# IDT (http://www.idtdna.com/calc/analyzer) uses
# (all biopython stuff from
# http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html )
# DNA/DNA 	                        +/- 1.4C (Allawi '97)
idtTable = NEAREST_N_ALLAWI_LUCIA_1997
# Divalent cation correction 	+/- 0.5C (Owczarzy '08) 
# Triphosphate correction 	        +/- 0.0C (Owczarzy '08)
# Note: monovalent must be done post; see GetIdtMeltingTemp, Below
idtSaltCorr = SALT_MONO_Owczarzy_2004
# default oligo concentration, nM
idtDnaOligoConc = 250
""" ibid:
"Oligo concentration [250nM] is assumed to be significantly larger
 (at least 6x) than concentration of the complementary target."
"""
idtOtherConc = idtDnaOligoConc/6

IdtParams = MeltingParams(Na=50, # 50mM by default
                          Tris=0,
                          dNTPs=0,
                          Mg = 0,
                          K=0,
                          saltcorr=idtSaltCorr,
                          nn_table=idtTable,
                          dnac1=idtDnaOligoConc,
                          dnac2=idtOtherConc,selfcomp=True)
    
# all possible salt corrrections
# from
#biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html#salt_correction
# note that some of these are redundant, so I just pick the ones that arent
# note that  5 through 7 are methods for something completely difference
SaltCorrections = [i for i in range(8)]
# all possible NN correction tables
# from
# http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-pysrc.html
NearestNeighborTablesDNA = [mt.DNA_NN1,mt.DNA_NN2,mt.DNA_NN3,
                            mt.DNA_NN4]

def MeltingTemperature(Sequence,
                       **kwargs):
    """
    Gets the melting temperature of a sequence. All concentrations in mM

    For saltcorr and nn_table (or just generally), see:
    http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html

    Args:
        Sequence: the sting we care about
       **kwargs: passed directly to melting temperature. See MeltingParams
    """
    mParams = MeltingParams(**kwargs)
    paramDict= mParams.__dict__
    return mt.Tm_NN(Seq(Sequence),**paramDict)

def GetAllMeltingTemperatures(sequence,**kwargs):
    """
    Gets the melting temperatures from *all possible* combinations of
    NN tables and salt corrections
    Args:
        sequence: sequence to use
        **kwargs: Arguments for MeltingTemperature, concentrations only
        (ie: dont try assing in saltcorr or nn_table, this does *all* of them)
    returns
        2-D matrix, element [i,j] is using salt correction i, nn_table j
        from the 'SaltCorrections' and 'NearestNeighborTablesDNA' tables
    """
    numNN = len(NearestNeighborTablesDNA)
    numSalt = len(SaltCorrections)
    toRet = np.zeros((numSalt,numNN))
    for i,saltcorr in enumerate(SaltCorrections):
        for j,nn_table in enumerate(NearestNeighborTablesDNA):
            toRet[i,j] = MeltingTemperature(sequence,saltcorr=saltcorr,
                                            nn_table=nn_table,**kwargs)
    return toRet
    

def GetIdtMeltingTemperature(sequence,**kwargs):
    """
    Gets the melting temperature of a primer, using what IDT does
    
    Args:
        sequence: see GetAllMeltingTemperatures
        **kwargs: see GetAllMeltingTemperatures
    Returns:
        the melting temperature, according to IDT
    """
    mParams = IdtParams
    mParamDict = mParams.__dict__
    rawMelt = MeltingTemperature(sequence,**mParamDict)
    # idt actually uses two different corrections, one for monovalent
    # See discussion at top, also:
    #biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html#Tm_NN
    # and bottom of http://www.idtdna.com/calc/analyzer
    corr = mt.salt_correction(method=SALT_DIV_Owczarzy_2008,
                              seq=sequence,**mParams.concdict())
    return 1/(1/(rawMelt)+corr)
