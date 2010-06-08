# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# My Stuff
from Entity import Entity


__doc__="The structure class, representing a macromolecular structure."


class Structure(Entity):
    """
    The Structure class contains a collection of Model instances.
    """
    def __init__(self, id):
        self.level="S"
        Entity.__init__(self, id)
    
    # Special methods
    
    def __repr__(self):
        return "<Structure id=%s>" % self.get_id()
    
    # Private methods
    
    def _sort(self, m1, m2):
        """Sort models.
        
        This sorting function sorts the Model instances in the Structure instance.
        The sorting is done based on the model id, which is a simple int that
        reflects the order of the models in the PDB file.
        
        Arguments:
        o m1, m2 - Model instances
        """
        return cmp(m1.get_id(), m2.get_id())
    
    # Public
    
    def get_chains(self):
        for m in self:
            for c in m:
                yield c
    
    def get_residues(self):
        for c in self.get_chains():
            for r in c:
                yield r
    
    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a

    def renumber_residues(self, start=1):
        """ Renumbers residues in a structure starting from 1. 
            Keeps numbering consistent with gaps.
        """
        
        for m in self:
            for c in m:
                fresidue_num = c.get_list()[0].get_id()[1]
                displace = start - fresidue_num
                for r in c:
                    r.id = (r.id[0], r.id[1]+displace, r.id[2])
    
    def search_ss_bonds(self, threshold=3):
        """ Searches S-S bonds based on distances between atoms in the structure (first model only).
            Average distance is 2.05A. Threshold is 3A default.
        """

        from itertools import combinations
        
        model = self.child_list[0]
        cysteines = filter( (lambda r: r.get_resname() == 'CYS'), model.get_residues() )
        
        pairs = combinations(cysteines, 2) # Iterator with pairs
        
        for cys_pair in pairs:
            if cys_pair[0]['SG'] - cys_pair[1]['SG'] < threshold:
                yield cys_pair
        
    def build_biological_unit(self):
        """
        """
        
        from copy import deepcopy # To copy structure object
        
        if self.header['biological_unit']:
            biomt_data = self.header['biological_unit']
        else:
            return "PDB File lacks appropriate REMARK 350 entries to build Biological Unit."
        
        for transformation in biomt_data:
            M = transformation[0] # Rotation Matrix
            T = transformation[1] # Translation Vector
                      
            for atom in rot_struct:
                atom.transform(M, T)
                