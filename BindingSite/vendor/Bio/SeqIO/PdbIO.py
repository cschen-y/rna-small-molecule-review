





"""Bio.SeqIO support for accessing sequences in PDB and mmCIF files."""

import collections
import warnings

from Bio import BiopythonParserWarning
from Bio.Data.PDBData import protein_letters_3to1
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator

_aa3to1_dict = {}
_aa3to1_dict.update(protein_letters_3to1)
_aa3to1_dict.update(protein_letters_3to1_extended)


def _res2aacode(residue, undef_code="X"):
    """Return the one-letter amino acid code from the residue name.

    Non-amino acid are returned as "X".
    """
    if isinstance(residue, str):
        return _aa3to1_dict.get(residue, undef_code)

    return _aa3to1_dict.get(residue.resname, undef_code)


def AtomIterator(pdb_id, structure):
    """Return SeqRecords from Structure objects.

    Base function for sequence parsers that read structures Bio.PDB parsers.

    Once a parser from Bio.PDB has been used to load a structure into a
    Bio.PDB.Structure.Structure object, there is no difference in how the
    sequence parser interprets the residue sequence. The functions in this
    module may be used by SeqIO modules wishing to parse sequences from lists
    of residues.

    Calling functions must pass a Bio.PDB.Structure.Structure object.


    See Bio.SeqIO.PdbIO.PdbAtomIterator and Bio.SeqIO.PdbIO.CifAtomIterator for
    details.
    """
    model = structure[0]
    for chn_id, chain in sorted(model.child_dict.items()):
        
        residues = [
            res
            for res in chain.get_unpacked_list()
            if _res2aacode(res.get_resname().upper()) != "X"
        ]
        if not residues:
            continue
        
        
        gaps = []
        rnumbers = [r.id[1] for r in residues]
        for i, rnum in enumerate(rnumbers[:-1]):
            if rnumbers[i + 1] != rnum + 1 and rnumbers[i + 1] != rnum:
                
                gaps.append((i + 1, rnum, rnumbers[i + 1]))
        if gaps:
            res_out = []
            prev_idx = 0
            for i, pregap, postgap in gaps:
                if postgap > pregap:
                    gapsize = postgap - pregap - 1
                    res_out.extend(_res2aacode(x) for x in residues[prev_idx:i])
                    prev_idx = i
                    res_out.append("X" * gapsize)
                else:
                    warnings.warn(
                        "Ignoring out-of-order residues after a gap",
                        BiopythonParserWarning,
                    )
                    
                    
                    res_out.extend(_res2aacode(x) for x in residues[prev_idx:i])
                    break
            else:
                
                res_out.extend(_res2aacode(x) for x in residues[prev_idx:])
        else:
            
            res_out = [_res2aacode(x) for x in residues]
        record_id = f"{pdb_id}:{chn_id}"
        
        
        
        

        record = SeqRecord(Seq("".join(res_out)), id=record_id, description=record_id)
        
        record.annotations["molecule_type"] = "protein"

        record.annotations["model"] = model.id
        record.annotations["chain"] = chain.id

        record.annotations["start"] = int(rnumbers[0])
        record.annotations["end"] = int(rnumbers[-1])
        yield record


class PdbSeqresIterator(SequenceIterator):
    """Parser for PDB files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over chains in a PDB file as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file

        The sequences are derived from the SEQRES lines in the
        PDB file header, not the atoms of the 3D structure.

        Specifically, these PDB records are handled: DBREF, DBREF1, DBREF2, SEQADV, SEQRES, MODRES

        See: http://www.wwpdb.org/documentation/format23/sect3.html

        This gets called internally via Bio.SeqIO for the SEQRES based interpretation
        of the PDB file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.pdb", "pdb-seqres"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...     print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Equivalently,

        >>> with open("PDB/1A8O.pdb") as handle:
        ...     for record in PdbSeqresIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...         print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Note the chain is recorded in the annotations dictionary, and any PDB DBREF
        lines are recorded in the database cross-references list.
        """
        super().__init__(source, fmt="PDB")
        self.cache = None

    def __next__(self):
        """Iterate over the records in the PDB file."""
        if self.cache is None:
            chains = collections.defaultdict(list)
            metadata = collections.defaultdict(list)

            rec_name = None
            for line in self.stream:
                rec_name = line[0:6].strip()
                if rec_name == "SEQRES":
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    chn_id = line[11]
                    
                    
                    residues = [_res2aacode(res) for res in line[19:].split()]
                    chains[chn_id].extend(residues)
                elif rec_name == "DBREF":
                    
                    pdb_id = line[7:11]
                    
                    chn_id = line[12]
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    database = line[26:32].strip()
                    
                    db_acc = line[33:41].strip()
                    
                    db_id_code = line[42:54].strip()
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    metadata[chn_id].append(
                        {
                            "pdb_id": pdb_id,
                            "database": database,
                            "db_acc": db_acc,
                            "db_id_code": db_id_code,
                        }
                    )
                elif rec_name == "DBREF1":
                    
                    pdb_id = line[7:11]
                    
                    chn_id = line[12]
                    
                    database = line[26:32].strip()
                    
                    db_id_code = line[47:67].strip()
                elif rec_name == "DBREF2":
                    
                    if pdb_id != line[7:11] or chn_id != line[12]:
                        raise ValueError("DBREF2 identifiers do not match")
                    
                    db_acc = line[18:40].strip()
                    metadata[chn_id].append(
                        {
                            "pdb_id": pdb_id,
                            "database": database,
                            "db_acc": db_acc,
                            "db_id_code": db_id_code,
                        }
                    )
                

            if rec_name is None:
                raise ValueError("Empty file.")

            self.cache = []
            for chn_id, residues in sorted(chains.items()):
                record = SeqRecord(Seq("".join(residues)))
                record.annotations = {"chain": chn_id}
                
                record.annotations["molecule_type"] = "protein"
                if chn_id in metadata:
                    m = metadata[chn_id][0]
                    record.id = record.name = f"{m['pdb_id']}:{chn_id}"
                    record.description = (
                        f"{m['database']}:{m['db_acc']} {m['db_id_code']}"
                    )
                    for melem in metadata[chn_id]:
                        record.dbxrefs.extend(
                            [
                                f"{melem['database']}:{melem['db_acc']}",
                                f"{melem['database']}:{melem['db_id_code']}",
                            ]
                        )
                else:
                    record.id = chn_id
                self.cache.append(record)
        try:
            record = self.cache.pop(0)
        except IndexError:
            raise StopIteration
        else:
            return record


class PdbAtomIterator(SequenceIterator):
    """Parser for structures in a PDB files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over structures in a PDB file as SeqRecord objects.

        Argument source is a file-like object or a path to a file.

        The sequences are derived from the 3D structure (ATOM records), not the
        SEQRES lines in the PDB file header.

        Unrecognised three letter amino acid codes (e.g. "CSD") from HETATM
        entries are converted to "X" in the sequence.

        In addition to information from the PDB header (which is the same for
        all records), the following chain specific information is placed in the
        annotation:

        record.annotations["residues"] = List of residue ID strings
        record.annotations["chain"] = Chain ID (typically A, B ,...)
        record.annotations["model"] = Model ID (typically zero)

        Where amino acids are missing from the structure, as indicated by
        residue numbering, the sequence is filled in with 'X' characters to
        match the size of the missing region, and  None is included as the
        corresponding entry in the list record.annotations["residues"].

        This function uses the Bio.PDB module to do most of the hard work. The
        annotation information could be improved but this extra parsing should
        be done in parse_pdb_header, not this module.

        This gets called internally via Bio.SeqIO for the atom based
        interpretation of the PDB file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.pdb", "pdb-atom"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        Equivalently,

        >>> with open("PDB/1A8O.pdb") as handle:
        ...     for record in PdbAtomIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        """
        

        
        from Bio.PDB.PDBParser import PDBParser

        structure = PDBParser().get_structure(None, source)
        pdb_id = structure.header["idcode"]
        if not pdb_id:
            warnings.warn(
                "'HEADER' line not found; can't determine PDB ID.",
                BiopythonParserWarning,
            )
            pdb_id = "????"

        records = []
        for record in AtomIterator(pdb_id, structure):
            
            record.annotations.update(structure.header)

            

            records.append(record)

        self.records = iter(records)

    def __next__(self):
        """Return the next SeqRecord from the PDB ATOM stream."""
        return next(self.records)


PDBX_POLY_SEQ_SCHEME_FIELDS = (
    "_pdbx_poly_seq_scheme.asym_id",  
    "_pdbx_poly_seq_scheme.mon_id",  
)

STRUCT_REF_FIELDS = (
    "_struct_ref.id",  
    "_struct_ref.db_name",  
    "_struct_ref.db_code",  
    "_struct_ref.pdbx_db_accession",  
)

STRUCT_REF_SEQ_FIELDS = (
    "_struct_ref_seq.ref_id",  
    "_struct_ref_seq.pdbx_PDB_id_code",  
    "_struct_ref_seq.pdbx_strand_id",  
)


class CifSeqresIterator(SequenceIterator):
    """Parser for chains in an mmCIF files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over chains in an mmCIF file as SeqRecord objects.

        Argument source is a file-like object or a path to a file.

        The sequences are derived from the _entity_poly_seq entries in the
        mmCIF file, not the atoms of the 3D structure.

        Specifically, these mmCIF records are handled: _pdbx_poly_seq_scheme
        and _struct_ref_seq. The _pdbx_poly_seq records contain sequence
        information, and the _struct_ref_seq records contain database
        cross-references.

        See:
        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/pdbx_poly_seq_scheme.html
        and
        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_ref_seq.html

        This gets called internally via Bio.SeqIO for the sequence-based
        interpretation of the mmCIF file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.cif", "cif-seqres"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...     print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Equivalently,

        >>> with open("PDB/1A8O.cif") as handle:
        ...     for record in CifSeqresIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...         print(record.dbxrefs)
        ...
        Record id 1A8O:A, chain A
        ['UNP:P12497', 'UNP:POL_HV1N5']

        Note the chain is recorded in the annotations dictionary, and any mmCIF
        _struct_ref_seq entries are recorded in the database cross-references
        list.
        """

        
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict

        chains = collections.defaultdict(list)
        metadata = collections.defaultdict(list)
        records = MMCIF2Dict(source)

        
        
        for field in (
            PDBX_POLY_SEQ_SCHEME_FIELDS + STRUCT_REF_SEQ_FIELDS + STRUCT_REF_FIELDS
        ):
            if field not in records:
                records[field] = []
            elif not isinstance(records[field], list):
                records[field] = [records[field]]

        for asym_id, mon_id in zip(
            records["_pdbx_poly_seq_scheme.asym_id"],
            records["_pdbx_poly_seq_scheme.mon_id"],
        ):
            mon_id_1l = _res2aacode(mon_id)
            chains[asym_id].append(mon_id_1l)

        
        struct_refs = {}
        for ref_id, db_name, db_code, db_acc in zip(
            records["_struct_ref.id"],
            records["_struct_ref.db_name"],
            records["_struct_ref.db_code"],
            records["_struct_ref.pdbx_db_accession"],
        ):
            struct_refs[ref_id] = {
                "database": db_name,
                "db_id_code": db_code,
                "db_acc": db_acc,
            }

        
        
        for ref_id, pdb_id, chain_id in zip(
            records["_struct_ref_seq.ref_id"],
            records["_struct_ref_seq.pdbx_PDB_id_code"],
            records["_struct_ref_seq.pdbx_strand_id"],
        ):
            struct_ref = struct_refs[ref_id]

            
            metadata[chain_id].append({"pdb_id": pdb_id})
            metadata[chain_id][-1].update(struct_ref)

        records = []  
        for chn_id, residues in sorted(chains.items()):
            record = SeqRecord(Seq("".join(residues)))
            record.annotations = {"chain": chn_id}
            
            record.annotations["molecule_type"] = "protein"
            if chn_id in metadata:
                m = metadata[chn_id][0]
                record.id = record.name = f"{m['pdb_id']}:{chn_id}"
                record.description = f"{m['database']}:{m['db_acc']} {m['db_id_code']}"
                for melem in metadata[chn_id]:
                    record.dbxrefs.extend(
                        [
                            f"{melem['database']}:{melem['db_acc']}",
                            f"{melem['database']}:{melem['db_id_code']}",
                        ]
                    )
            else:
                record.id = chn_id
            records.append(record)  

        self.records = iter(records)

    def __next__(self):
        """Return the next SeqRecord from the mmCIF chain stream."""
        return next(self.records)


class CifAtomIterator(SequenceIterator):
    """Parser for structures in an mmCIF files."""

    modes = "t"

    def __init__(self, source: _TextIOSource) -> None:
        """Iterate over structures in an mmCIF file as SeqRecord objects.

        Argument source is a file-like object or a path to a file.

        The sequences are derived from the 3D structure (_atom_site.* fields)
        in the mmCIF file.

        Unrecognised three letter amino acid codes (e.g. "CSD") from HETATM
        entries are converted to "X" in the sequence.

        In addition to information from the PDB header (which is the same for
        all records), the following chain specific information is placed in the
        annotation:

        record.annotations["residues"] = List of residue ID strings
        record.annotations["chain"] = Chain ID (typically A, B ,...)
        record.annotations["model"] = Model ID (typically zero)

        Where amino acids are missing from the structure, as indicated by
        residue numbering, the sequence is filled in with 'X' characters to
        match the size of the missing region, and  None is included as the
        corresponding entry in the list record.annotations["residues"].

        This function uses the Bio.PDB module to do most of the hard work. The
        annotation information could be improved but this extra parsing should
        be done in parse_pdb_header, not this module.

        This gets called internally via Bio.SeqIO for the atom based
        interpretation of the PDB file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("PDB/1A8O.cif", "cif-atom"):
        ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        Equivalently,

        >>> with open("PDB/1A8O.cif") as handle:
        ...     for record in CifAtomIterator(handle):
        ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        ...
        Record id 1A8O:A, chain A

        """
        

        
        from Bio.PDB.MMCIFParser import MMCIFParser

        structure = MMCIFParser().get_structure(None, source)
        pdb_id = structure.header["idcode"]
        if not pdb_id:
            warnings.warn("Could not determine the PDB ID.", BiopythonParserWarning)
            pdb_id = "????"
        self.records = AtomIterator(pdb_id, structure)

    def __next__(self):
        """Return the next SeqRecord from the mmCIF ATOM stream."""
        return next(self.records)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
