# Copyright 2010 by Joao Rodrigues.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading MODELLER PIR format files as
# SeqRecord objects.  The code is based on Bio.SeqIO.PirIO.

"""Bio.SeqIO support for MODELLER's PIR file format.

You are expected to use this module via the Bio.SeqIO functions, or if
the file contains a sequence alignment, optionally via Bio.AlignIO instead.

This format is used by MODELLER to describe both sequences and alignments.

The file format is described online at:
http://www.salilab.org/modeller/9v8/manual/node454.html

An example file in this format would be:

>P1;5fd1
structureX:5fd1:1    :A:106  :A:ferredoxin:Azotobacter vinelandii: 1.90: 0.19
AFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA
EVWPNITEKKDPLPDAEDWDGVKGKLQHLER*

>P1;1fdx
sequence:1fdx:1    : :54   : :ferredoxin:Peptococcus aerogenes: 2.00:-1.00
AYVINDSC--IACGACKPECPVNIIQGS--IYAIDADSCIDCGSCASVCPVGAPNPED-----------------
-------------------------------*


As with the FASTA format, each record starts with a line begining with ">"
character.  There is then a two letter sequence type, P1, a semi colon, and 
the protein identification code. This code corresponds to the sequence code.

The second line contains information necessary to extract atomic coordinates
from the original PDB coordinate set. It is composed of 10 fields separated by 
colon characters (;). The fields are as follows:

1 - 
A specification of whether or not 3D structure is available and of the type 
of the method used to obtain the structure: structureX, X-ray; structureN, NMR; 
structureM, model; sequence, sequence. Only structure is also a valid value.
    
2 - 
The PDB filename or code. Can be a path. 

3-6 - 
The residue and chain identifiers (see below) for the first (fields 3-4) and 
last residue (fields 5-6) of the sequence in the subsequent lines. 

A residue number is unspecified when a blank character or a dot, `.', is given.
A chain id is not specified when a dot, `.', is given.

7 - 
Protein name. Optional.

8 - 
Source of the protein. Optional.

9 - 
Resolution of the crystallographic analysis. Optional.

10 - 
R-factor of the crystallographic analysis. Optional.


The remaining lines contain the sequence itself,
terminating in an asterisk.

The alignment file can contain blank lines between the protein entries.

Comment lines can occur outside protein entries. 
They must begin with the identifiers `C;' or `R;'.

"""

from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#This is a generator function!
def mPirIterator(handle):
    """Generator function to iterate over mPIR records (as SeqRecord objects).

    Arguments:
        handle - input file
    
    """
    
    allowed_specifications = {  'structureX': 'X-Ray Structure', 
                                'structureN': 'NMR Structure', 
                                'structureM': 'Model', 
                                'sequence': 'Sequence', 
                                'structure': 'Generic Structure'}
    
    allowed_chain_ids = set('abcdefghijklmnopqrstuvwxyz0123456789@')
    
    #Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return #Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        
        # First line: protein code
        if line[0:4] != ">P1;":
            raise ValueError(\
                "Records in mPIR files should start with '>P1;'")
        
        protein_code = line[4:].strip()
        
        # Second line: extraction info
        # Fields separated by ':' character.
        
        line = handle.readline().strip()
        fields = line.split(':')
        
        if len(fields) != 10:
            raise ValueError(\
                "Second line of each record must have 10 fields\
                 separated by the colon (:) character. Fields may be empty.")
        
        if fields[0] not in allowed_specifications:
            raise ValueError(\
                "Unrecognized specification (%s)" %fields[0])

        record_type = allowed_specifications[fields[0]]
        pdb_code = fields[1]
        
        # Residue and Chain Boundaries
        init_res = fields[2].strip()
        init_chain = fields[3].strip().lower()
        end_res = fields[4].strip()
        end_chain = fields[5].strip().lower()
        
        if not ( (init_res in ["", "FIRST"] or init_res.isdigit())
                 and 
                 (end_res in ["", "END", "LAST", "+nn"] or end_res.isdigit()) ):
                
            raise ValueError(\
                "Unrecognized residue identifier: Initial (%s), End (%s)" 
                %(init_res, end_res))
        
        if not ((init_chain in allowed_chain_ids or init_chain == "")
                and
                (end_chain in allowed_chain_ids or end_chain == "")):

            raise ValueError(\
                "Unrecognized chain identifiers: Initial (%s), End (%s)" 
                %(init_chain, end_chain))
        
        # Optional Fields
        
        protein_name = fields[6]
        source_organism = fields[7]
        resolution = fields[8]
        r_factor = fields[9]
        
        # Until *, sequence.
            
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            #Remove trailing whitespace, and any internal spaces
            lines.append(line.rstrip().replace(" ",""))
            line = handle.readline()
        seq = "".join(lines)
        if seq[-1] != "*":
            #Note the * terminator is present on nucleotide sequences too,
            #it is not a stop codon!
            raise ValueError(\
                "Sequences in mPIR files should include a * terminator!")
            
        #Return the record and then continue...
        record = SeqRecord(Seq(seq[:-1], generic_protein),
                           id = protein_code, name = protein_code,
                           description = protein_name)
        record.annotations["record_type"] = record_type
        record.annotations["initial_residue"] = init_res
        record.annotations["end_residue"] = end_res
        record.annotations["initial_chain"] = init_chain
        record.annotations["end_chain"] = end_chain
        record.annotations["source_organism"] = source_organism
        record.annotations["resolution"] = resolution
        record.annotations["r_factor"] = r_factor
        
        yield record
  
        if not line : return #StopIteration
    assert False, "Should not reach this line"