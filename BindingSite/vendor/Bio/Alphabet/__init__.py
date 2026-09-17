







"""Alphabets were previously used to declare sequence type and letters (OBSOLETE).

The design of Bio.Alphabet included a number of historic design choices
which, with the benefit of hindsight, were regretable. Bio.Alphabet was
therefore removed from Biopython in release 1.78. Instead, the molecule type is
included as an annotation on SeqRecords where appropriate.

Please see https://biopython.org/wiki/Alphabet for examples showing how to
transition from Bio.Alphabet to molecule type annotations.
"""

raise ImportError(
    "Bio.Alphabet has been removed from Biopython. In many cases, the alphabet can simply be ignored and removed from scripts. In a few cases, you may need to specify the ``molecule_type`` as an annotation on a SeqRecord for your script to work correctly. Please see https://biopython.org/wiki/Alphabet for more information."
)
