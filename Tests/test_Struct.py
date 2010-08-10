# Copyright 2010 by Joao Rodrigues.  All rights reserved.
# 
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Struct module."""

import os
import unittest
import warnings

try:
    from numpy.random import random
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "Install NumPy if you want to use Bio.PDB.")

from Bio import Struct
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning

class test_Protein(unittest.TestCase):
    
    def setUp(self):
        warnings.resetwarnings()
        warnings.simplefilter('ignore', PDBConstructionWarning)
        
        structure = Struct.read('PDB/1A8O.pdb')
        self.p = structure.as_protein()
        
    def test_Hydrogenate(self):
        """Hydrogenation Algorithm """
        
        from Bio.Struct import Hydrogenate as H
        
        p = self.p
        
        protonate = H.Hydrogenate_Protein()
        protonate.add_hydrogens(p)
        
        self.assertEqual(1037, len( [i for i in p.get_atoms() ]))
        
    
    def test_SSSearch(self):
        """Sulfide Bridge Search """
        
        p = self.p
        
        ss_bonds = p.search_ss_bonds()
        
        self.assertEqual(1, len( [i for i in ss_bonds ]))
    
    def test_CheckMissing(self):
        """Check Missing Atoms """
        
        p = self.p
        missing = p.check_missing_atoms()
        self.assertEqual(0, len( missing.keys() ))
    
    def test_CoarseGrained(self):
        """Coarse Grained Models """
        
        p = self.p
        
        # CA
        cg_model = p.coarse_grain()
        self.assertEqual(66, len([i for i in cg_model.get_atoms()]))

        # ENCAD
        cg_model = p.coarse_grain('ENCAD_3P')
        self.assertEqual(194, len([i for i in cg_model.get_atoms()]))


        # MARTINI
        cg_model = p.coarse_grain('MARTINI')
        self.assertEqual(142, len([i for i in cg_model.get_atoms()]))
        
        
if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)