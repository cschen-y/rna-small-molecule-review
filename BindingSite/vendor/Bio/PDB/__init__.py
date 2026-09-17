






"""Classes that deal with macromolecular crystal structures.

Includes: PDB and mmCIF parsers, a Structure class, a module to keep a local
copy of the PDB up-to-date, selective IO of PDB files, etc.

Original Author: Thomas Hamelryck.
Contributions by:
- Peter Cock
- Joe Greener
- Rob Miller
- Lenna X. Peterson
- Joao Rodrigues
- Kristian Rother
- Eric Talevich
- and many others.
"""

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install NumPy if you want to use Bio.PDB. See http://www.numpy.org/"
    ) from None




from . import Selection


from .cealign import CEAligner


from .Dice import extract



from .DSSP import DSSP
from .DSSP import make_dssp_dict


from .FragmentMapper import FragmentMapper
from .HSExposure import ExposureCN


from .HSExposure import HSExposureCA
from .HSExposure import HSExposureCB
from .mmcifio import MMCIFIO
from .MMCIFParser import FastMMCIFParser
from .MMCIFParser import MMCIFParser


from .parse_pdb_header import parse_pdb_header


from .PDBIO import PDBIO
from .PDBIO import Select


from .PDBList import PDBList
from .PDBMLParser import PDBMLParser
from .PDBParser import PDBParser
from .Polypeptide import CaPPBuilder
from .Polypeptide import is_aa
from .Polypeptide import is_nucleic


from .Polypeptide import PPBuilder
from .Polypeptide import standard_aa_names
from .ResidueDepth import get_surface



from .ResidueDepth import ResidueDepth


from .StructureAlignment import StructureAlignment


from .Superimposer import Superimposer
from .vectors import calc_angle
from .vectors import calc_dihedral
from .vectors import m2rotaxis
from .vectors import refmat
from .vectors import rotaxis
from .vectors import rotaxis2m
from .vectors import rotmat


from .vectors import Vector
from .vectors import vector_to_axis



try:
    from .NeighborSearch import NeighborSearch
except ImportError:
    pass



try:
    from .SASA import ShrakeRupley
except ImportError:
    pass
