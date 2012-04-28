# Copyright (C) 2011, Joao Rodrigues (j.rodrigues@uu.nl)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrapper and for structural file formats.

This API follows the same semantics as Biopython's SeqIO and AlignIO.
"""

# For 'with' on Python/Jython 2.5
from __future__ import with_statement
__docformat__ = "restructuredtext en"

import os.path
import warnings

from Bio import File

# Write-related imports
from Bio.PDB import PDBIO

supported_o_formats = {
        'pdb': PDBIO
        }

# Read-related imports
from Bio.PDB import PDBParser

supported_i_formats = {
        'pdb': PDBParser
        }

# Check for an mmCIF parser. It may be missing, either because of an unbuilt C
# extension (standard lex/yacc) or the pure-Python fallback which depends on the
# third-party package PLY.
try:
    from Bio.PDB import MMCIFParser
except ImportError:
    # mmCIF format isn't available, but don't bother the user about it until
    # it's requested as a format in the read() function (below).
    pass
else:
    supported_i_formats['mmcif'] = MMCIFParser


# Helpers

def _file_to_id(path_or_handle):
    """Take file name without the extension."""
    if hasattr(path_or_handle, 'name'):
        # It's a 'file' object that knows its own path
        full_name = path_or_handle.name
    elif isinstance(path_or_handle, basestring):
        # It's a path to a file
        full_name = path_or_handle
    else:
        # File-like stream, e.g. StringIO, pipe or network handle
        return "<unknown structure>"
    # Full path -- remove filename extension and directory path
    no_ext = os.path.splitext(full_name)[0]
    return os.path.basename(no_ext)


# Public Functions

def read(infile, format='pdb', id=None, **kwargs):
    """Parse a molecular structure from the given file or handle.

    Wrapper for Bio.PDB parsers.

    Other parser-specific arguments can be passed as keyworded arguments.

    :Parameters:
        infile : string or file-like object
            Path to file, or file-like handle
        format : string (optional)
            File format. One of 'pdb' or 'mmcif'; default 'pdb'.
        id : string (optional)
            Structure Id. Typically the PDB ID, but can be anything.
    """
    if format.lower() == 'mmcif' and 'mmcif' not in supported_i_formats:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError(
            "Module MMCIParser is missing because of unmet dependencies "
            "(unbuilt C extension or missing PLY package).")

    if id is None:
        id = _file_to_id(infile)
    # Get format and instantiate
    Parser = supported_i_formats[format.lower()](**kwargs)
    get_structure = getattr(Parser, 'get_structure')
    return get_structure(id, infile)


def write(structure, outfile, format='pdb', **kwargs):
    """Write a structure to a file or handle in the given format.

    Other writer-specific arguments can be passed as keyworded arguments.

    :Parameters:
        structure : Bio.PDB.Structure instance
            Structure Data (SMCRA representation)
        outfile : string or file-like object
            Path to file, or file-like handle
        format : string (optional)
            File format. Currently, only 'pdb' is supported.
    """
    w = supported_o_formats[format.lower()]
    Writer = w(**kwargs)

    with File.as_handle(outfile, 'w') as handle:
        set_structure = getattr(Writer, 'set_structure')
        set_structure(structure)
        save = getattr(Writer, 'save')
        return save(handle)

