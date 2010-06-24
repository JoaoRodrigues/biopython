# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This module is heavily based on PyMol's code and also in MMTK's code.
# Similarities are not a coincidence.
# PyMol: chempy/protein.py chempy/place.py
# MMTK: Proteins.py / ChemicalObjects.py
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Vector Operations

import numpy as N

# BioPython modules
from Bio.PDB.Entity import Entity
from Bio.PDB.Atom import Atom

# Bond information (Taken from PyMol)
import protein_residues
import protein_amber
import bond_amber

class Hydrogenate_Protein:

    def __init__(self, u_input):
        
        # Load Structure & Perform routines
        self.nh_structure = u_input
        self.exclude_ss_bonded_cysteines()
        
        # FF and Tmplt pre-load
        
        self.tmpl = {   'cter': protein_residues.c_terminal,
                        'nter': protein_residues.n_terminal,
                        'nrml': protein_residues.normal
                    }
        self.ffld = {   'cter': protein_amber.c_terminal,
                        'nter': protein_amber.n_terminal,
                        'nrml': protein_amber.normal
                    }

        # Constants (Taken from PyMol)
        self.N_TERMINAL_ATOMS = set(['HT','HT1','HT2','HT3','H1','H2','H3',
                                  '1H','2H','3H','1HT','2HT','3HT'])

        self.C_TERMINAL_ATOMS = set(['OXT','O2','OT1','OT2'])

        # Geometries adapted from MMTK (add Length and Angles)
        self.geometries =   {   'C':    {   3: 'Tetrahedral',
                                            2: 'Trigonal Planar',
                                            1: 'Linear'
                                        },
                                'N':    {   3: 'Tetrahedral',
                                            2: 'Trigonal Planar',
                                            1: 'Linear'
                                        },
                                'O':    {   1: 'Linear'},
                                'S':    {   1: 'Linear'},
                            }
        
    
    def exclude_ss_bonded_cysteines(self):
        
        # Pre-compute ss bonds to discard cysteines later
        ss_bonds =  self.nh_structure.search_ss_bonds()
        self.ss_cysteines = []
        for cys_pair in ss_bonds:
            self.ss_cysteines += list(cys_pair)
        
        return self.ss_cysteines
        
    def find_missing_hydrogens(self):
        """
        Takes as input an Entity-based class (S,M,C,R) or a list of Residues.
        Iterates over the object to find missing atoms in protein residues.
        Uses protein_amber.py, bond_amber.py, and protein_residues.py to find those missing.
        Adds them to the structure with coordinates x=y=z=666 and serial number 0.
        """

        if isinstance(self.nh_structure, Entity): # Structure, Model, Chain, Residues
            if self.nh_structure.level in ['S', 'M', 'C']:
                residue_list = list(self.nh_structure.get_residues())
            else:
                residue_list = [self.nh_structure]
        elif hasattr(self.nh_structure, '__iter__'): # List of Residues?
            if not [i for i in self.nh_structure if not i.level == 'R']: # Every single object is a residue
                residue_list = self.nh_structure
        else: # Some other weirdo object
            raise ValueError('Can only search for missing atoms in Entity-based classes (S,M,C,R).')

            
        self.missing_per_residue = {}

        for residue in residue_list:
            
            self.missing_per_residue[residue] = {}
            
            if residue in self.ss_cysteines:
                exclude_sh = True
            else:
                exclude_sh = False
            
            atom_list = residue.child_dict
            names = set(atom_list.keys())

            # Set Template & Force Field
            if names.intersection(self.C_TERMINAL_ATOMS):
                selection = 'cter'
            elif names.intersection(self.N_TERMINAL_ATOMS):
                selection = 'nter'
            else:
                selection = 'nrml'

            # Define stuff
            tmpl = self.tmpl[selection]
            ffld = self.ffld[selection]

            # Go!
            if not tmpl.has_key(residue.resname):
                raise ValueError("Unknown Residue Type: %s" %residue.resname)
            else:
                bonds = tmpl[residue.resname]['bonds']

                for pair in bonds.keys(): # Returns tuples of bonded atoms e.g. (N, CA)
                    a, b = pair
                    if a in names and b not in names: # b is missing
                        if a == 'SG' and exclude_sh:
                            continue
                        if not ( b[0] == 'H' or ( b[1] == 'H' and b[0].isdigit()) ): # Hx and dHx
                            continue
                        
                        # Update missing
                        if self.missing_per_residue[residue].has_key(atom_list[a]):
                            self.missing_per_residue[residue][atom_list[a]].append(b)
                        else:
                            self.missing_per_residue[residue][atom_list[a]] = [b]
                        
                    elif b in names and a not in names: # a is missing
                        if b == 'SG' and exclude_sh:
                            continue
                        if not ( a[0] == 'H' or ( a[1] == 'H' and a[0].isdigit()) ): # Hx and dHx
                            continue

                        if self.missing_per_residue[residue].has_key(atom_list[b]):
                            self.missing_per_residue[residue][atom_list[b]].append(a)
                        else:
                            self.missing_per_residue[residue][atom_list[b]] = [a]
                         
                    else: # Both exist
                        pass

        return 
            
    def find_anchor(self, atom, tmpl):
        """
        Given a heavy atom, finds another heavy atom it is connected to.
        Necessary for H positioning.
        Returns Atom object.
        """
        residue = atom.parent
        bonds = tmpl[residue.resname]['bonds']
        
        for pair in bonds:
            a,b = pair
            if ( a[0] == 'H' or ( a[0].isdigit() and a[1] == 'H') ) or ( b[0] == 'H' or ( b[0].isdigit() and b[1] == 'H') ):
                pass # Ignore Hs
            elif a == atom.name:
                anchor = [i for i in residue.child_list if i.name == b][0]
                return anchor
            elif b == atom.name:
                anchor = [i for i in residue.child_list if i.name == a][0]
                return anchor
        
        return 0   
                
        
    def build_hydrogens(self, bond_field=bond_amber):
        """
        Places missing hydrogens atoms based on simple geometric criteria.
        """
             
        for residue in self.missing_per_residue:
                        
            print "Building hydrogens on residue:", residue.resname
            
            need = [ [], [], [], [] ]
            
            

            
            
            
            
            
            
            
            
            # # Iterate..
            # for heavy_atom in ha_list:
            #     
            #     h_atoms = ha_list[heavy_atom]
            #     
            #     nH = len(h_atoms) # Number of Hs to add
            #     element = heavy_atom.element # Element of Core Atom (C,N,S,O)
            #     
            #     if not self.geometries.has_key(element):
            #         try:
            #             element = ffld[(residue.resname, heavy_atom.name)]['type'][0] # Try to recover element from ffld
            #         except KeyError:
            #             raise ValueError('Residue or Atom not recognized: (%s, %s)' %(residue.resname, heavy_atom.name))
            #     
            #     geometry = self.geometries[element][nH] # Tetrahedral, etc...
            # 
            #         
            #     anchor =  self.find_anchor(heavy_atom, tmpl) # Second heavy atom
            # 
            #     print "Heavy Atom: %s / Element is %s / Geometry is %s / Anchor is %s" %(heavy_atom, element, geometry, anchor)
            #     print nH
            #     try:
            #         new_atom = geometry(heavy_atom, anchor, None)
            #         yield new_atom
            #     except:
            #         print 'No can do'
            #     
            #     # 
                    
                
                
                
        
        