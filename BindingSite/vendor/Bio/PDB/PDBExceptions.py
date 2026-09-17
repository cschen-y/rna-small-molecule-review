






"""Some Bio.PDB-specific exceptions."""

from Bio import BiopythonWarning



class PDBException(Exception):
    """Define class PDBException."""




class PDBConstructionException(Exception):
    """Define class PDBConstructionException."""


class PDBConstructionWarning(BiopythonWarning):
    """Define class PDBConstructionWarning."""



class PDBIOException(Exception):
    """Define class PDBIOException."""


class PDBIOWarning(BiopythonWarning):
    """Define class PDBIOWarning."""
