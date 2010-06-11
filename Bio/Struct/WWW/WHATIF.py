# Copyright 2010 by Joao Rodrigues.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with the WWW version of WHATIF
servers provided at http://swift.cmbi.kun.nl/.

Functions:
add_hydrogens        Add Hydrogens to a protein using the HTOPO service.
"""

import urllib
import os, tempfile

# Utility functions
def smcra_to_pdb(smcra_object, temp_dir='/tmp/'):
    """
    Since all servers work with PDB files, we can't use the SMCRA structure directly.
    Using PDBIO & tempfile to make it work seamlessly.
    """
    
    temp_path = tempfile.mktemp( '.pdb', dir=temp_dir )
    
    io = PDBIO()
    io.set_structure(smcra_object)
    io.save(temp_path)
    
    f = open(temp_path, 'r')
    pdb_data = f.read()
    f.close()
    
    os.remove(temp_path)
    
    return pdb_data

def pdb_to_smcra(pdb_contents, temp_dir='/tmp/'):
    
    temp_path = tempfile.mktemp( '.pdb', dir=temp_dir )
    f = open(temp_path, 'r')
    f.write(pdb_contents)
    f.close()
        
    p = PDBParser()
    smcra_object = p.get_structure(temp_path)
 
    os.remove(temp_path)
    
    return smcra_object

# Dummy REST webservice to Test WHATIF servers

def whatif_test(): # Check if service is up
    pass

# HTOPO Service to Add Hydrogens

def add_hydrogens(data):
    """ From the web server: http://swift.cmbi.ru.nl/servers/html/htopo.html
        All missing protons will be added to the structure.
    """
    
    if not isinstance(data, str): # Hope it's a SMCRA object..
        data = smcra_to_pdb(data)
        smcra_data = 1
    else:
        smcra_data = 0
    
    h_added = whatif(data, "htopo")
    
    if smcra_data:
        return pdb_to_smcra(h_added)
    else:
        return h_added

# Waiting for G.Vriend to add REST access