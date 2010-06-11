# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB.Structure import Structure
from Bio.PDB import Entity

class Protein(Structure):
    
    def __init__(self, protein):

        Structure.__init__(self, protein.id)
        
        self.full_id = protein.full_id
        self.child_list = protein.child_list
        self.child_dict = protein.child_dict
        self.parent = protein.parent
        self.xtra = protein.xtra
        
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