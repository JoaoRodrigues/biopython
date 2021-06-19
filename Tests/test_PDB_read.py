# Copyright 2021 by Simon Duerr.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for PDB module."""

import gzip
import sys
import unittest
import warnings

from pathlib import Path

from Bio import PDB

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class Test_PDB_read(unittest.TestCase):
    """Tests for Bio.PDB.read function."""

    def test_get_file_extension(self):
        """Test private get_file_extension function."""
        fn_list = ["1a.pdb", "1a.pdb.gz", "1a.renamed.pdb.gz", "1a.renamed.pdb"]
        fn_list += [Path(fn) for fn in fn_list]  # test Path objects too.
        for fn in fn_list:
            self.assertEqual(PDB._get_file_extension(fn), "pdb")

    def test_PDB_read_pathlike(self):
        """Test PDB.read() from pathlib.Path."""
        fp = Path("PDB/1A8O.pdb")
        fp_gz = Path("PDB/1A8O.pdb.gz")

        _ = PDB.read(fp)
        _ = PDB.read(fp, fmt="pdb")
        _ = PDB.read(fp_gz)
        _ = PDB.read(fp_gz, fmt="pdb")

    def test_PDB_read_str(self):
        """Test PDB.read() from string."""
        fn = "PDB/1A8O.pdb"
        fn_gz = "PDB/1A8O.pdb.gz"

        _ = PDB.read(fn)
        _ = PDB.read(fn, fmt="pdb")
        _ = PDB.read(fn_gz)
        _ = PDB.read(fn_gz, fmt="pdb")

    def test_PDB_read_handle(self):
        """Test PDB.read() from open file handle."""

        with open("PDB/1A8O.pdb") as fp:
            _ = PDB.read(fp, fmt="pdb")

        with gzip.open("PDB/1A8O.pdb.gz", mode="rt") as fp:
            _ = PDB.read(fp, fmt="pdb")

        with self.assertRaises(TypeError):
            with open("PDB/1A8O.pdb") as fp:
                _ = PDB.read(fp)

    # MMTF Format
    def test_PDB_read_mmtf(self):
        """Test PDB.read() with MMTF files."""

        try:
            from Bio.PDB.mmtf import MMTFParser
        except ImportError:
            self.skipTest("dependency 'mmtf' is not installed")
        else:
            _ = PDB.read("PDB/1A8O.mmtf")

    def test_PDB_read_mmtf_handle_error(self):
        """Test PDB.read() raises error when reading MMTF from handle."""

        try:
            from Bio.PDB.mmtf import MMTFParser
        except ImportError:
            self.skipTest("dependency 'mmtf' is not installed")
        else:
            # MMTF reader does not accept handles
            with self.assertRaises(TypeError):
                with open("PDB/1A8O.mmtf") as handle:
                    _ = PDB.read(handle, format="mmtf")

    def test_PDB_read_kwargs(self):
        """Test PDB.read() with parser kwargs."""

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            _ = PDB.read("PDB/occupancy.pdb", PERMISSIVE=True)
            self.assertEqual(len(w), 3)

        # With quiet mode no warnings should be printed
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            _ = PDB.read("PDB/occupancy.pdb", QUIET=True, PERMISSIVE=True)
            self.assertEqual(len(w), 0)

    def test_PDB_read_format_error(self):
        """Test PDB.read() fmt errors."""

        # uppercase format should not raise an error
        _ = PDB.read("PDB/1A8O.pdb", fmt="PDB")

        # empty format
        with self.assertRaises(ValueError):
            _ = structure = PDB.read("PDB/1A8O.pdb", fmt="")
        # unknown format
        with self.assertRaises(ValueError):
            _ = structure = PDB.read("PDB/1A8O.pdb", fmt="blabla")

    def test_PDB_read_extract_id(self):
        """Test PDB.read() extracting 'id' from structure header."""

        # No header, no ID
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = PDB.read("PDB/disordered.pdb")
            self.assertEqual(structure.id, None)

        # ID is read from PDB header
        structure = PDB.read("PDB/1A8O.pdb")
        self.assertEqual(structure.id, "1A8O")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
