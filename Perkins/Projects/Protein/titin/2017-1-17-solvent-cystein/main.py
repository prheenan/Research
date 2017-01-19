import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored

import numpy as np
# god forgive me for absolute paths
sys.path.append("/Users/patrickheenan/src_prh/")
from GeneralUtil.python import CheckpointUtilities


def GetArea(identifier,solvent_area,dot_density=4,solvent_radius=1.4):
    """
    Args:
        identifier: the idenfier for what we are getting the area of;
        passed to cmd.get_area

        dot_density: the density of dots, between 0 and 1
        
        solvent_area,: if True, get the solvent-accessible area. Otherwise,
        just get the surface area
    
        solvent_radius: radius of the solvent, in Angstroms (def is water)
    """
    cmd.set('dot_solvent', int(solvent_area))
    cmd.set('dot_density', dot_density)
    cmd.set('solvent_radius', solvent_radius)
    return cmd.get_area(identifier)

def GetAllAreas(name_cystines):
    # loop through, get the areas and accessible areas
    sasa_per_residue,area_per_residue = [],[]
    for i in stored.residues:
        solvent_area = GetArea("resi %s" % i,solvent_area=True)
        molecular_area = GetArea("resi %s" % i,solvent_area=False)
        sasa_per_residue.append(solvent_area)
        area_per_residue.append(molecular_area)
    return sasa_per_residue,area_per_residue


def run():
    # reset the session, make the background less awful
    cmd.delete("all")
    cmd.bg_color('grey80')
    base="/Users/patrickheenan/src_prh/Research/Perkins/Projects/Protein/" +\
        "titin/2017-1-17-solvent-cystein/"
    cmd.load(base +'1waa.pdb')  # use the name of your pdb file
    maximum_accessible_area_angstroms_squared = 3
    stored.residues = []
    # get all the cystines, only really care about them
    name_cystines = "cys"
    cmd.select(name_cystines,"resn " + name_cystines)
    cmd.iterate(name_cystines, 'stored.residues.append(resi)')
    ids = [float(r) for r in stored.residues]
    sasa_per_residue,area_per_residue= \
        CheckpointUtilities.getCheckpoint(base + "out.pkl",
                                          GetAllAreas,False,name_cystines)
    accessible_tmp = np.where(sasa_per_residue > 
                              maximum_accessible_area_angstroms_squared)[0]
    ids_greater_than = [stored.residues[i] for i in accessible_tmp]
    # make all the cystines spheres
    cmd.show("spheres",name_cystines)
    accessible_name = "accessible_name"
    # color the accessible ones black
    cmd.select(accessible_name,"resi " + "+".join(ids_greater_than))
    cmd.color("black",accessible_name)
    # get the names of all the accessible ones...
    stored.accessible_names, stored.accessible_chains = [],[]
    cmd.iterate(name_cystines, 'stored.accessible_names.append(resn)')
    cmd.iterate(name_cystines, 'stored.accessible_chains.append(chain)')
    fraction =np.array(sasa_per_residue)/np.array(area_per_residue)
    # save out the data...
    header_data = [ ["chain",stored.accessible_chains],
                    ["residue name" ,stored.accessible_names],
                    ["residue number",ids],
                    ["accessible area (AA^2)",sasa_per_residue],
                    ["total residue area (AA^2)",area_per_residue],
                    ["fraction accessible",fraction]]
    header = ",".join([h[0] for h in header_data])
    data = np.array([h[1] for h in header_data],dtype=np.object).T
    print(data)
    np.savetxt(fname=base + "out.csv",X=data,fmt="%s",delimiter=",",
               header=header)


run()

