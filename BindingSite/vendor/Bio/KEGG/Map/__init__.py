







"""Load KEGG Pathway maps for use with the Biopython Pathway module.

The pathway maps are in the format::

    RXXXXX:[X.X.X.X:] A + 2 B <=> C
    RXXXXX:[X.X.X.X:] 3C <=> 2 D + E
    ...

where RXXXXX is a five-digit reaction id, and X.X.X.X is the optional
EC number of the enzyme that catalyze the reaction.
"""

from Bio.Pathway import Reaction


def parse(handle):
    """Parse a KEGG pathway map."""
    for line in handle:
        data, catalysts, reaction = line.split(":")
        catalysts = [(catalysts,)]
        reactants = {}
        before, after = reaction.split("<=>")
        compounds = before.split(" + ")
        for compound in compounds:
            compound = compound.strip()
            try:
                number, compound = compound.split()
                number = -int(number)
            except ValueError:
                number = -1
            reactants[compound] = number
        compounds = after.split(" + ")
        for compound in compounds:
            compound = compound.strip()
            try:
                number, compound = compound.split()
                number = int(number)
            except ValueError:
                number = +1
            reactants[compound] = number
        yield Reaction(reactants, catalysts, True, data)
