# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This module is heavily based on PyMol's code.
# Similarities are not a coincidence.
# PyMol: chempy/protein.py chempy/place.py
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Vector Operations

from cpv import *

# BioPython modules
from Bio.PDB.Entity import Entity
from Bio.PDB.Atom import Atom

# Bond information (Taken from PyMol)
import protein_residues
import protein_amber
import bond_amber

TET_TAN = 1.41
TRI_TAN = 1.732

class Hydrogenate_Protein:

    def __init__(self, u_input, forcefield=protein_amber, template=protein_residues):
        
        # Load Structure & Perform routines
        self.nh_structure = u_input
        self._exclude_ss_bonded_cysteines()
        
        # FF and Tmplt pre-load
        
        self.tmpl = {   'cter': template.c_terminal,
                        'nter': template.n_terminal,
                        'nrml': template.normal
                    }
        self.ffld = {   'cter': forcefield.c_terminal,
                        'nter': forcefield.n_terminal,
                        'nrml': forcefield.normal
                    }

        # Constants (Taken from PyMol)
        self.N_TERMINAL_ATOMS = set(['HT','HT1','HT2','HT3','H1','H2','H3',
                                  '1H','2H','3H','1HT','2HT','3HT'])

        self.C_TERMINAL_ATOMS = set(['OXT','O2','OT1','OT2'])
        
        # Protonation Methods
        
        self.protonation_methods = {
                                    1: self._add_1,
                                    2: self._add_2,
                                    3: self._add_3,
                                    4: self._add_4
                                    }
    
    def _build_bonding_network(self):
        """
        Evaluates atoms per residue for missing and known bonded partners.
        Based on bond_amber.
        A better alternative would be to iterate over the entire list of residues and 
        use NeighborSearch to probe neighbors for atom X in residue i, i-1 and i+1
        """
        
        self.bonds = {} # [Residue][Atom]: [ [missing], [bonded] ]
        self.selection = {} # [Residue]: 'nrml'/'nter'/'cter'
        
        for residue in self.nh_structure.get_residues():
            
            bond_dict = self.bonds[residue] = {}
            
            atom_dict = residue.child_dict
            atom_names = set(atom_dict.keys())
            
            # Pre-Populate Dictionary
            
            for name in atom_names:
                bond_dict[name] = [ [], [] ]
            
            # Define Template
            if atom_names.intersection(self.C_TERMINAL_ATOMS):
                selection = 'cter'
            elif atom_names.intersection(self.N_TERMINAL_ATOMS):
                selection = 'nter'
            else:
                selection = 'nrml'
            
            tmpl = self.tmpl[selection]
            self.selection[residue] = selection # For place_hs
            
            # Iterate Template Bonds and record info
            
            if not tmpl.has_key(residue.resname):
                raise ValueError("Unknown Residue Type: %s" %residue.resname)
            
            template_bonds = tmpl[residue.resname]['bonds']
            
            for bond in template_bonds.keys():
                
                a1, a2 = bond
                
                if a1 in atom_names and not a2 in atom_names:
                    bond_dict[a1][0].append(a2)
                elif a1 not in atom_names and a2 in atom_names:
                    bond_dict[a2][0].append(a1)
                else: # 
                    bond_dict[a1][1].append(atom_dict[a2])
                    bond_dict[a2][1].append(atom_dict[a1])   
    
    def _exclude_ss_bonded_cysteines(self):
        
        # Pre-compute ss bonds to discard cystines for H-adding.
        ss_bonds =  self.nh_structure.search_ss_bonds()
        for cys_pair in ss_bonds:
            cys1, cys2 = cys_pair
            
            cys1.resname = 'CYX'
            cys2.resname = 'CYX'
    
    def _find_secondary_anchors(self, residue, heavy_atom, anchor):
        """
        Searches through the bond network for atoms bound to the anchor.
        Returns a secondary and tertiary anchors.
        Example, for CA, returns C and O.
        """
        
        for secondary in self.bonds[residue][anchor.name][1]:
            for tertiary in self.bonds[residue][secondary.name][1]:
                if tertiary.name != heavy_atom.name and tertiary.name != anchor.name:
                    return (secondary,tertiary)
        
        return None
        
    def place_hydrogens(self, bondfield=bond_amber):
        
        # define bondfield
        self.bondfield = bondfield
        
        last_count = -1
        
        # Determine how many Hs each HA needs
        
        for residue in sorted(self.bonds.keys(), key=lambda x: x.get_id()[1]):
            incomplete_atoms = self.bonds[residue]
            print
            print residue.resname
            for atom in incomplete_atoms:
                missing_atoms = incomplete_atoms[atom][0]
                if len(missing_atoms):                
                
                    #near = self._find_secondary_anchors(residue, atom, anchor)
                    #print near.name
                    #h_add = self.protonation_methods[len(missing_atoms)]
                    miss = missing_atoms[0]
                    print atom
                    if len(missing_atoms) == 1:
                        h_coord = self._add_1(miss, residue.child_dict[atom], incomplete_atoms[atom][1])
                        new_at = Atom(name=miss, coord=h_coord, bfactor=1.0, occupancy=0.0, 
                                        altloc=' ', fullname=miss, serial_number=random.randint(5000,9999), element='H' )
                        residue.add(new_at)
                        
                    # elif len(missing_atoms) == 2:
                    #     coordinates = self._add_2(missing_atoms, residue.child_dict[atom], incomplete_atoms[atom][1])
                    #     for name in coordinates:
                    #         new_at = Atom(name=name, coord=coordinates[name], bfactor=1.0, occupancy=0.0, 
                    #                         altloc=' ', fullname=name, serial_number=random.randint(5000,9999), element='H')
                    #         # residue.add(new_at)
                    #         print new_at.name,
                    #         print new_at.coord
                            
    def _add_1(self, hydrogen_name, heavy_atom, bonds):
        """
        Adds one proton to an atom.
        """
        
        residue = heavy_atom.parent
        ffld = self.ffld[self.selection[residue]]
        bnd_len = self.bondfield.length
        anchor = bonds[0]
                
        # If not linear
        if self.bondfield.nonlinear.has_key(ffld[(residue.resname, heavy_atom.name)]['type']):
            bonded = self._find_secondary_anchors(residue, heavy_atom, anchor) # Tuple of two atoms

            if bonded:
                if self.bondfield.planer.has_key(ffld[(residue.resname, anchor.name)]['type']): # Phenolic hydrogens, etc.
                    secondary_anchor = bonded[0]
                   
                    p0 = heavy_atom.coord - anchor.coord
                    d2 = secondary_anchor.coord - anchor.coord
                    p1 = normalize(cross_product(d2,p0))
                    p2 = normalize(cross_product(p0,p1))                     
                    v = scale(p2,TRI_TAN)
                    v = normalize(add(p0,v))
                    
                    hydrogen_coord = add(heavy_atom.coord,
                                        scale(v, 
                                              bnd_len[(ffld[(residue.resname, heavy_atom.name)]['type'], 
                                              ffld[(residue.resname, hydrogen_name)]['type'])]    ))
                
                else: # Ser, Cys, Thr hydroxyl hydrogens
                    secondary_anchor = bonded[0]
                    v = anchor.coord - secondary_anchor.coord
                    hydrogen_coord = add(heavy_atom.coord,
                                        scale(v, 
                                              bnd_len[(ffld[(residue.resname, heavy_atom.name)]['type'], 
                                              ffld[(residue.resname, hydrogen_name)]['type'])]    ))
                    
            elif len(bonds):
                d2 = [1.0,0,0]
                p0 = heavy_atom.coord -bonds[0].coord
                p1 = normalize(cross_product(d2,p0))
                v = scale(p1,TET_TAN)
                v = normalize(add(p0,v))
                
                hydrogen_coord = add(heavy_atom.coord,
                                     scale(  v, 
                                             bnd_len[(ffld[(residue.resname, heavy_atom.name)]['type'], 
                                             ffld[(residue.resname, hydrogen_name)]['type'])]    ))
            else:
                hydrogen_coord = random_sphere(heavy_atom.coord,
                                               bnd_len[ (ffld[(residue.resname, heavy_atom.name)]['type'],
                                                         ffld[(residue.resname, hydrogen_name)]['type']) ])

        elif len(bonds): # linear sum...amide, tbu, etc
            v = [0.0,0.0,0.0]
            if heavy_atom.name == 'N': # Fix to get previous atom O from peptide bond. Ugly.
                prev_res = list(residue.get_id())
                prev_res[1] -= 1
                prev_res = tuple(prev_res)
                if residue.parent.child_dict.has_key(prev_res):
                    prev_res = residue.parent.child_dict[prev_res]
                    bonds.append(prev_res.child_dict['O'])
            for b in bonds:
                d = heavy_atom.coord - b.coord
                v = add(v, d)
            v = normalize(v)
            hydrogen_coord = add(heavy_atom.coord,
                                 scale(v,
                                 bnd_len[(ffld[(residue.resname, heavy_atom.name)]['type'],
                                          ffld[(residue.resname, hydrogen_name)]['type']) ]))
        else:
            hydrogen_coord = random_sphere(heavy_atom.coord,
                                           bnd_len[ (ffld[(residue.resname, heavy_atom.name)]['type'],
                                                     ffld[(residue.resname, hydrogen_name)]['type']) ])
        
        return hydrogen_coord
        
    def _add_2(self, miss, atom, anchor, bonded):
        
        coord = {} # Returns two coordinate sets
        
        residue = atom.parent
        sel = self.selection[residue]
        at1 = atom
        at2 = miss[0]
        know = bonded
        bnd_len = self.bondfield.length

        if self.bondfield.planer.has_key(self.ffld[sel][(residue.resname, at1.name)]['type']): # guanido, etc
            near = self._find_secondary_anchors(residue, at1, anchor)
            if near: # 1-4 present
                print "Planar Near"
                at3 = anchor
                at4 = near
                d1 = sub(at1.coord,at3.coord)
                p0 = normalize(d1)
                d2 = sub(at4.coord,at3.coord)
                p1 = normalize(cross_product(d2,p0))
                p2 = normalize(cross_product(p0,p1))
                v = scale(p2,TRI_TAN)
                v = normalize(add(p0,v))
                coord[at2] = add(at1.coord,scale(v,
                  bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))                                                         
                at2 = miss[1]
                v = scale(p2,-TRI_TAN)
                v = normalize(add(p0,v))
                coord[at2] = add(at1.coord,scale(v,
                  bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))
            
            elif len(know): # no 1-4 found
                print "Planar not/near"
                d2 = [1.0,0,0]
                at3 = anchor
                d1 = sub(at1.coord,at3.coord)
                p0 = normalize(d1)                  
                p1 = normalize(cross_product(d2,p0))
                p2 = normalize(cross_product(p0,p1))
                v = scale(p2,TRI_TAN)
                v = normalize(add(p0,v))
                coord[at2] = add(at1.coord,scale(v,
                  bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))
                at2 = miss[1]
                v = scale(p2,-TRI_TAN)
                v = normalize(add(p0,v))
                coord[at2] = add(at1.coord,scale(v,
                  bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))
            else:
                coord = random_sphere(at1.coord,
                    bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])])
        
        elif len(know)>=2: # simple tetrahedral
            print "Simple T"
            at3 = anchor
            at4 = bonded[1]
            v = [0.0,0.0,0.0]
            d1 = sub(at1.coord,at3.coord)
            d2 = sub(at1.coord,at4.coord)
            v = add(normalize(d1),normalize(d2))
            p0 = normalize(v)
            p1 = normalize(cross_product(d2,p0))
            v = scale(p1,TET_TAN)
            v = normalize(add(p0,v))
            coord[at2] = add(at1.coord,scale(v,
                    bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))
            at2 = miss[1]               
            v = scale(p1,-TET_TAN)
            v = normalize(add(p0,v))
            coord[at2] = add(at1.coord,scale(v,
                    bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))
        else:
            if len(know): # sulfonamide? 
                print "Sulfonamide"
                d2 = [1.0,0,0]
                at3 = anchor
                d1 = sub(at1.coord,at3.coord)
                p0 = normalize(d1)                                    
                p1 = normalize(cross_product(d2,p0))
                v = scale(p1,TET_TAN)
                v = normalize(add(p0,v))
                coord[at2] = add(at1.coord,scale(v,
                    bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))
            else: # blind
                coord[at2] = random_sphere(at1.coord,
                    bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])])
            at4=at2
            at2=miss[1]
            v = [0.0,0.0,0.0]
            d1 = sub(at1.coord,at3.coord)
            d2 = sub(at1.coord,at4.coord)
            v = add(normalize(d1),normalize(d2))
            p0 = normalize(v)
            p1 = normalize(cross_product(d2,p0))
            v = scale(p1,TET_TAN)
            v = normalize(add(p0,v))
            coord[at2] = add(at1.coord,scale(v,
                bnd_len[(self.ffld[sel][(residue.resname, at1.name)]['type'],self.ffld[sel][(residue.resname, at2)]['type'])]))
        return coord
        
    def _add_3(self):
        pass
    def _add_4(self):
        pass