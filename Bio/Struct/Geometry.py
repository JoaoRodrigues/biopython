# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB import Entity

def center_of_mass(entity, geometric=False):
    """
    Returns gravitic or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    Defaults to Gravitic.
    """
    
    if isinstance(entity, Entity.Entity): # Structure, Model, Chain, Residue
        atom_list = entity.get_atoms()
    elif hasattr(entity, '__iter__') and filter(lambda x: x.level == 'A', entity): # List of Atoms
        atom_list = entity
    else: # Some other weirdo object
        raise ValueError('Center of Mass can only be calculated from the following objects:\n \
                            Structure, Model, Chain, Residue, list of Atoms.')
    
    masses = []
    pos = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]
    
    for atom in atom_list:
        masses.append(atom.mass)
        
        for i, coord in enumerate(atom.coord.tolist()):
            pos[i].append(coord)
    
    if '?' in set(masses): # If there is a single atom with undefined mass, default to Geometric
        geometric = True
        import warnings
        from PDBExceptions import PDBConstructionWarning
        warnings.warn("Some atoms have undefined masses. Calculating Geometric Center of Mass instead.",
                      PDBConstructionWarning)
    
    if geometric:
        return map( lambda y: sum(y)/len(masses), pos)
    else:
        w_pos = map(lambda x: [x[i]*masses[i] for i in range(len(x))], pos)
        return map(lambda y: sum(y)/sum(masses), w_pos)    