# Copyright (C) 2002,2020 Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes that deal with macromolecular crystal structures.

Includes: PDB and mmCIF parsers, a Structure class, a module to keep a local
copy of the PDB up-to-date, selective IO of PDB files, etc.

Original Author: Thomas Hamelryck.
Contributions by:
- Peter Cock
- Simon Duerr
- Joe Greener
- Rob Miller
- Lenna X. Peterson
- Joao Rodrigues
- Kristian Rother
- Eric Talevich
- and many others.
"""

import gzip
from pathlib import Path, PurePath

# Get a Structure object from a PDB file using available parsers
from .PDBParser import PDBParser
from .MMCIFParser import MMCIFParser
from .MMCIFParser import FastMMCIFParser

try:
    from .mmtf import MMTFParser
    _have_mmtf = True
except ImportError:
    _have_mmtf = False

# Download from the PDB
from .PDBList import PDBList

# Parse PDB header directly
from .parse_pdb_header import parse_pdb_header

# Find connected polypeptides in a Structure
from .Polypeptide import PPBuilder, CaPPBuilder, is_aa, standard_aa_names

# This is also useful :-)
from Bio.Data.SCOPData import protein_letters_3to1

# IO of PDB files (including flexible selective output)
from .PDBIO import PDBIO, Select
from .mmcifio import MMCIFIO

# Some methods to eg. get a list of Residues
# from a list of Atoms.
from . import Selection

# Superimpose atom sets
from .Superimposer import Superimposer

# 3D vector class
from .vectors import Vector, calc_angle, calc_dihedral, refmat, rotmat, rotaxis
from .vectors import vector_to_axis, m2rotaxis, rotaxis2m

# Alignment module
from .StructureAlignment import StructureAlignment

# DSSP handle
# (secondary structure and solvent accessible area calculation)
from .DSSP import DSSP, make_dssp_dict

# Residue depth:
# distance of residue atoms from solvent accessible surface
from .ResidueDepth import ResidueDepth, get_surface

# Calculation of Half Sphere Solvent Exposure
from .HSExposure import HSExposureCA, HSExposureCB, ExposureCN

# Kolodny et al.'s backbone libraries
from .FragmentMapper import FragmentMapper

# Write out chain(start-end) to PDB file
from .Dice import extract

# Fast atom neighbor search
# Depends on kdtrees C module
try:
    from .NeighborSearch import NeighborSearch
except ImportError:
    pass

# Native Shrake-Rupley algorithm for SASA calculations.
# Depends on kdtrees C module
try:
    from .SASA import ShrakeRupley
except ImportError:
    pass

_fmt_to_parser = {
    # PDB
    "pdb": PDBParser,
    "pdb.gz": PDBParser,
    "pdb1": PDBParser,
    "pdb1.gz": PDBParser,
    "ent": PDBParser,
    "ent.gz": PDBParser,
    # mmCIF
    "mmcif": MMCIFParser,
    "cif": MMCIFParser,
    "cif.gz": MMCIFParser,
}

if _have_mmtf:
    _fmt_to_parser["mmtf"] = MMTFParser
    _fmt_to_parser["mmtf.gz"] = MMTFParser


def _get_file_extension(filename):
    """Returns the file extension, ignoring compression suffixes (such as .gz).
    """
    exts = [ext[1:].lower() for ext in Path(filename).suffixes]
    if exts[-1] == 'gz':
        return exts[-2]
    return exts[-1]


def _as_handle_gzip(handle):
    """Open str, bytes or Path like object - gzipped or not."""
    try:
        with open(handle, "rb") as fp:
            magic = fp.read(2)
            fp.seek(0)
            if magic == b"\x1f\x8b":
                return gzip.open(handle, mode="rt")
            else:
                return open(handle)
    except TypeError:
        # should already be handle
        return handle


def read(handle, fmt=None, **kwargs):
    """Create a Structure from an open file or a file name.


    Arguments:
     - handle    - handle to the open file, or file name as a string.
     - fmt       - string defining the file format. If None (default), will try
                   to guess the file format from the file extension.

    Examples:

    You can read a structure from a file name:

    >>> from Bio import PDB
    >>> filename = "../Tests/PDB/1A8O.pdb"
    >>> structure = PDB.read(filename)
    >>> structure = PDB.read(filename, fmt="pdb")

    or from a pathlib Path object:

    >>> import pathlib
    >>> filepath = pathlib.Path("../Tests/PDB/1A8O.pdb")
    >>> structure = PDB.read(filepath)

    or from an open file handle by specifying the format explicitly:

    >>> from Bio import PDB
    >>> with open("../Tests/PDB/1A8O.pdb.gz") as handle:
    ...    structure = PDB.read(handle, fmt="pdb")

    You can also pass any optional arguments to read() as you would to the
    individual parsers:

    >>> from Bio import PDB
    >>> structure = PDB.read("../Tests/PDB/1A8O.pdb", QUIET=True)
    """
    # Validate format or try guessing from file name.
    if isinstance(handle, (str, PurePath)) and fmt is None:
        fmt = _get_file_extension(handle)
    elif fmt is None:
        raise TypeError(
            "Argument 'fmt' is required when reading from an open file handle."
        )

    fmt = fmt.lower()
    if fmt not in _fmt_to_parser:
        raise ValueError(
            f"'{fmt}' is not one of the following supported formats: "
            f"{', '.join(_fmt_to_parser.keys())}"
        )

    # Instantiate parser
    parser = _fmt_to_parser[fmt](**kwargs)

    # MMTF needs a separate parser, since its get_structure() takes only one
    # positional argument and does not support file handles.
    if _have_mmtf and isinstance(parser, MMTFParser):
        if not isinstance(handle, str):
            raise TypeError(
                "MMTFParser does not support reading from an open file handle."
            )
        parser = MMTFParser()
        return parser.get_structure(handle)

    # PDB and MMCIF parsers
    with _as_handle_gzip(handle) as fp:
        structure = parser.get_structure(None, fp)
        # If the structure has header info and header contains an idcode
        # try setting it to that value.
        if hasattr(structure, "header"):
            structure.id = structure.header.get("idcode").strip() or None

    return structure
