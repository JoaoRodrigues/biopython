# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from copy import deepcopy

from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
from Bio.PDB.Polypeptide import is_aa

class Protein(Structure):
    
    def __init__(self, id):
        
        Structure.__init__(self, id)
        
    @classmethod
    def from_structure(cls, original, filter_residues):
        
        P = cls(original.id)
        P.full_id = original.full_id
        
        for child in original.child_dict.values():
            copycat = deepcopy(child)
            P.add(copycat)
        
        # Discriminate non-residues (is_aa function)
        if filter_residues:
            for model in P:
                for chain in model:
                    map(chain.detach_child, [res.id for res in chain.child_list if not is_aa(res)])
            for chain in P.get_chains(): # Remove empty chains
                if not len(chain.child_list):
                    model.detach_child(chain.id)
                    
        P.header = deepcopy(original.header)
        P.xtra = deepcopy(original.xtra)
        
        return P
          
    def search_ss_bonds(self, threshold=3.0):
        """ Searches S-S bonds based on distances between atoms in the structure (first model only).
            Average distance is 2.05A. Threshold is 3A default.
            Returns iterator with tuples of residues.
        """

        from itertools import combinations
        
        model = self.child_list[0]
        cysteines = filter( (lambda r: r.get_resname() == 'CYS'), model.get_residues() )

        pairs = combinations(cysteines, 2) # Iterator with pairs

        for cys_pair in pairs:
            if cys_pair[0]['SG'] - cys_pair[1]['SG'] < threshold:
                yield cys_pair
    
def coarse_grain(ori_entity, cg_type="CA"):
    """ 
        Reduces the protein structure complexity to a few atoms per residue.
        Arguments:
        - ori_entity: Structure, Model, Chain, or List of Residues.
        Parameters:
        - cg_type:      CA (Ca-only) [Default]
                        3pt (Ca + O + Side-chain C.o.M)
        Returns a new structure object.
    """
    
    from Bio.Struct.Geometry import center_of_mass
    
    entity = deepcopy(ori_entity)
    
    if isinstance(entity, Entity): # S,M,C
        residue_list = entity.get_residues()
    elif hasattr(entity, "__iter__") and filter(lambda x: x.level == 'R', entity): # list of R
        residue_list = entity
    else: # Some other weirdo object
        raise ValueError('You can only use coarse grain on Structure, Model, Chain, or list of Residue(s).')
        
    if cg_type == "CA":
        for residue in residue_list:
            for atom_name in residue.child_dict.keys():
                if atom_name != 'CA':
                    residue.detach_child(atom_name)
    
    elif cg_type == "3pt":
        
        from Bio.PDB.Atom import Atom # To add C.O.M
        
        for residue in residue_list:

            side_chain_atoms = []
            for atom in residue.child_dict.keys():
                if atom in ['N', 'C']:
                    residue.detach_child(atom)
                elif atom not in ['CA', 'O']:
                    side_chain_atoms.append(residue[atom])
                    residue.detach_child(atom)
            if len(side_chain_atoms) == 1: # Alanine
                sc_com = side_chain_atoms[0]
                sc_com.name = 'CMA'
                residue.add(sc_com)                    
            elif len(side_chain_atoms) > 1:
                sc_com_coord = center_of_mass(side_chain_atoms)
                if residue.resname in ['GLU', 'ASP']:
                    sc_com = Atom('CMA', sc_com_coord, 0, 0, ' ', 'CMA', 0, 'O')
                elif residue.resname in ['LYS', 'ARG', 'HIS']:
                    sc_com = Atom('CMA', sc_com_coord, 0, 0, ' ', 'CMA', 0, 'N')
                elif residue.resname == 'CYS':
                    sc_com = Atom('CMA', sc_com_coord, 0, 0, ' ', 'CMA', 0, 'S')
                else:
                    sc_com = Atom('CMA', sc_com_coord, 0, 0, ' ', 'CMA', 0, 'C')

                residue.add(sc_com)
                          
    else:
        raise ValueError("Invalid complexity value: %s\nMust either be 1 or 3." %complexity)
        
    return entity  