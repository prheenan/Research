# (c) patrick raymond heenan, 2017, all right reserved
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored
import numpy as np,pickle

def GetArea(identifier,solvent_area,dot_density=4,solvent_radius=1.4):
    """
    Args:
        identifier: the idenfier for what we are getting the area of;
        passed to cmd.get_area

        dot_density: the density of dots, between 0 and 4
        
        solvent_area,: if True, get the solvent-accessible area. Otherwise,
        just get the surface area
    
        solvent_radius: radius of the solvent, in Angstroms (def is water)
    """
    cmd.set('dot_solvent', int(solvent_area))
    cmd.set('dot_density', dot_density)
    cmd.set('solvent_radius', solvent_radius)
    return cmd.get_area(identifier)

def GetAllAreas(chains,residues,names):
    """
    gets the area of each individual atom in the residue and chain

    Args:
        chain: which chain to use
        residues: numbers for the residues
        names: name of the atom within a given residue
    Returns:
        surface area and solvent area for the given residues 
    """
    # loop through, get the areas and accessible areas
    sasa_per_residue,area_per_residue = [],[]
    for chain,residue,name in zip(chains,residues,names):
        my_str = "chain {:s} & resi {:s} & name {:s}".format(chain,residue,name)
        solvent_area = GetArea(my_str,solvent_area=True)
        molecular_area = GetArea(my_str,solvent_area=False)
        sasa_per_residue.append(solvent_area)
        area_per_residue.append(molecular_area)
    return sasa_per_residue,area_per_residue


def run():
    """
    Gets all the accessible areas of the cystines
    
    TO RUN:
    
    (1) open pymol
    (2) in the command window, type:
    run /Users/patrickheenan/src_prh/Research/Perkins/Projects/Protein/titin/2017-1-17-solvent-cystein/main.py
    (replacing whereever the file is. Note you must also change 'base' below
    to point to where thre file lives. 
    (3) Then this will spit out a CSV file with accessibility information on 
    cystines. it also colors the cystines red/green if they are/arent accessible
    (in pymol)
    """
    # reset the session, make the background less awful
    cmd.delete("all")
    cmd.bg_color('grey80')
    # need to replace this with wherever the pdb file is...
    # must be an abolute path
    base="/Users/patrickheenan/src_prh/Research/Perkins/Projects/Protein/" +\
        "titin/2017-1-17-solvent-cystein/"
    """
    Note: lyle sent me chain A (I think) of titin IG27
    """
    cmd.load(base +'1waa.pdb')
    maximum_accessible_area_angstroms_squared = 3
    stored.residues,stored.names,stored.chains = [],[],[]
    # get all the cystines, only really care about them
    name_cystines = "cys"
    cmd.select(name_cystines,"resn " + name_cystines)
    cmd.iterate(name_cystines, 'stored.chains.append(chain)')
    cmd.iterate(name_cystines, 'stored.residues.append(resi)')
    cmd.iterate(name_cystines, 'stored.names.append(name)')
    ids = [float(r) for r in stored.residues]
    cache_file = base + "cache.pkl"
    if (os.path.isfile(cache_file)):
        with open(cache_file,'rb') as f:
            sasa_per_residue,area_per_residue = pickle.load(f)
    else:
        sasa_per_residue,area_per_residue= GetAllAreas(stored.chains,
                                                       stored.residues,
                                                       stored.names)
        with open(cache_file,'wb') as f:
            pickle.dump((sasa_per_residue,area_per_residue),f,
                        protocol=pickle.HIGHEST_PROTOCOL)
    accessible_tmp = np.where(sasa_per_residue > 
                              maximum_accessible_area_angstroms_squared)[0]
    ids_greater_than = [stored.residues[i] for i in accessible_tmp]
    # make all the cystines spheres
    cmd.show("spheres",name_cystines)
    cmd.color("green",name_cystines)
    accessible_name = "accessible_cys"
    # color the accessible ones red
    cmd.select(accessible_name,"resi " + "+".join(ids_greater_than))
    cmd.color("red",accessible_name)
    thiol_name = "thiol"
    cmd.select(thiol_name, "resn cys & name SG")
    cmd.show("spheres","thiol")
    cmd.color("yellow","thiol")
    # get the names of all the accessible ones...
    stored.cys_residue_names, stored.cys_atomic_names,stored.cys_chains = \
        [],[],[]
    cmd.iterate(name_cystines, 'stored.cys_residue_names.append(resn)')
    cmd.iterate(name_cystines, 'stored.cys_chains.append(chain)')
    cmd.iterate(name_cystines, 'stored.cys_atomic_names.append(name)')
    """
    The names property in pymol is a
    "list of up to 4-letter codes for atoms in proteins or nucleic acids"
    (see http://pymol.sourceforge.net/newman/user/S0220commands.html)

    The code is available also here: 
    http://chemistry.umeche.maine.edu/Modeling/NAMD8.html

    more details here: 
    PROTEIN DATA BANK, 
    ATOMIC COORDINATE AND BIBLIOGRAPHIC ENTRY FORMAT DESCRIPTION, 1992
    cdn.rcsb.org/rcsb-pdb//file_formats/pdb/pdbguide2.2/PDB_format_1992.pdf

    (1) greek letters are what you would think
    (2) Cystine goes like:
    
    C_alpha -- C_beta -- S_gamma

    or, in PDB notation:

    CA -- CB -- SG
    (see Table III, a of:
    Abbreviations and Symbols for the Description of the Conformation of 
    Polypeptide Chains. European Journal of Biochemistry 17, 193-201 (1970).)

    too bad, looks like the thiol (S_gamma) is accessible)
    """
    fraction =np.array(sasa_per_residue)/np.array(area_per_residue)
    # save out the data...
    header_data = [ ["chain",stored.cys_chains],
                    ["residue name" ,stored.cys_residue_names],
                    ["residue number",ids],
                    ["atom name",stored.cys_atomic_names],
                    ["accessible area (AA^2)",sasa_per_residue],
                    ["total residue area (AA^2)",area_per_residue],
                    ["fraction accessible",fraction]]
    header = ",".join([h[0] for h in header_data])
    data = np.array([h[1] for h in header_data],dtype=np.object).T
    np.savetxt(fname=base + "out.csv",X=data,fmt="%s",delimiter=",",
               header=header)


run()

