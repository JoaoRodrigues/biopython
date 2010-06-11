# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for dealing with PDB structures.
"""

import os

def read( handle, id=None ):

    from Bio.PDB import PDBParser
    
    if not id:
        id = os.path.basename(handle).split('.')[0] # Get from filename
    
    p = PDBParser()
    s = p.get_structure(id, handle)

    return s