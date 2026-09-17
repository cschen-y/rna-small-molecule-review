






r"""Bio.SearchIO support for Bill Pearson's FASTA tools.

This module adds support for parsing FASTA outputs. FASTA is a suite of
programs that finds regions of local or global similarity between protein
or nucleotide sequences, either by searching databases or identifying
local duplications.

Bio.SearchIO.FastaIO was tested on the following FASTA flavors and versions:

    - flavors: fasta, ssearch, tfastx
    - versions: 35, 36

Other flavors and/or versions may introduce some bugs. Please file a bug report
if you see such problems to Biopython's bug tracker.

More information on FASTA are available through these links:

    - Website: http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml
    - User guide: http://fasta.bioch.virginia.edu/fasta_www2/fasta_guide.pdf


Supported Formats
=================

Bio.SearchIO.FastaIO supports parsing and indexing FASTA outputs triggered by
the -m 10 flag. Other formats that mimic other programs (e.g. the BLAST tabular
format using the -m 8 flag) may be parseable but using SearchIO's other parsers
(in this case, using the 'blast-tab' parser).


fasta-m10
=========

Note that in FASTA -m 10 outputs, HSPs from different strands are considered to
be from different hits. They are listed as two separate entries in the hit
table. FastaIO recognizes this and will group HSPs with the same hit ID into a
single Hit object, regardless of strand.

FASTA also sometimes output extra sequences adjacent to the HSP match. These
extra sequences are discarded by FastaIO. Only regions containing the actual
sequence match are extracted.

The following object attributes are provided:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | program                 | FASTA flavor                     |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | full length of query sequence    |
|                 +-------------------------+----------------------------------+
|                 | target                  | target search database           |
|                 +-------------------------+----------------------------------+
|                 | version                 | FASTA version                    |
+-----------------+-------------------------+----------------------------------+
| Hit             | seq_len                 | full length of the hit sequence  |
+-----------------+-------------------------+----------------------------------+
| HSP             | bitscore                | \*_bits line                     |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | \*_expect line                   |
|                 +-------------------------+----------------------------------+
|                 | ident_pct               | \*_ident line                    |
|                 +-------------------------+----------------------------------+
|                 | init1_score             | \*_init1 line                    |
|                 +-------------------------+----------------------------------+
|                 | initn_score             | \*_initn line                    |
|                 +-------------------------+----------------------------------+
|                 | opt_score               | \*_opt line, \*_s-w opt line     |
|                 +-------------------------+----------------------------------+
|                 | pos_pct                 | \*_sim line                      |
|                 +-------------------------+----------------------------------+
|                 | sw_score                | \*_score line                    |
|                 +-------------------------+----------------------------------+
|                 | z_score                 | \*_z-score line                  |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | aln_annotation          | al_cons block, if present        |
| (also via HSP)  +-------------------------+----------------------------------+
|                 | hit                     | hit sequence                     |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | hit sequence end coordinate      |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit sequence start coordinate    |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query                   | query sequence                   |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query sequence end coordinate    |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query sequence start coordinate  |
|                 +-------------------------+----------------------------------+
|                 | query_strand            | query sequence strand            |
+-----------------+-------------------------+----------------------------------+

"""

import re

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import Hit
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult

__all__ = ("FastaM10Parser", "FastaM10Indexer")




_RE_FLAVS = re.compile(r"t?fast[afmsxy]|pr[sf][sx]|lalign|[gs]?[glso]search")

_PTR_ID_DESC_SEQLEN = r">>>(.+?)\s+(.*?) *- (\d+) (?:aa|nt)\s*$"
_RE_ID_DESC_SEQLEN = re.compile(_PTR_ID_DESC_SEQLEN)
_RE_ID_DESC_SEQLEN_IDX = re.compile(_PTR_ID_DESC_SEQLEN.encode())

_RE_ATTR = re.compile(r"^; [a-z]+(_[ \w-]+):\s+(.*)$")

_RE_START_EXC = re.compile(r"^-*")
_RE_END_EXC = re.compile(r"-*$")


_HSP_ATTR_MAP = {
    "_initn": ("initn_score", int),
    "_init1": ("init1_score", int),
    "_opt": ("opt_score", int),
    "_s-w opt": ("opt_score", int),
    "_z-score": ("z_score", float),
    "_bits": ("bitscore", float),
    "_expect": ("evalue", float),
    "_score": ("sw_score", int),
    "_ident": ("ident_pct", float),
    "_sim": ("pos_pct", float),
}


_STATE_NONE = 0
_STATE_QUERY_BLOCK = 1
_STATE_HIT_BLOCK = 2
_STATE_CONS_BLOCK = 3


def _set_qresult_hits(qresult, hit_rows=()):
    """Append Hits without alignments into QueryResults (PRIVATE)."""
    for hit_row in hit_rows:
        hit_id, remainder = hit_row.split(" ", 1)
        
        
        
        
        
        
        if hit_id not in qresult:
            frag = HSPFragment(hit_id, qresult.id)
            hsp = HSP([frag])
            hit = Hit([hsp])
            qresult.append(hit)

    return qresult


def _set_hsp_seqs(hsp, parsed, program):
    """Set HSPs sequences (PRIVATE).

    :param hsp: HSP whose properties will be set
    :type hsp: HSP
    :param parsed: parsed values of the HSP attributes
    :type parsed: dictionary {string: object}
    :param program: program name
    :type program: string

    """
    
    start = 0
    for seq_type in ("hit", "query"):
        if "tfast" not in program:
            pseq = parsed[seq_type]
            
            
            start, stop = _get_aln_slice_coords(pseq)
            start_adj = len(re.search(_RE_START_EXC, pseq["seq"]).group(0))
            stop_adj = len(re.search(_RE_END_EXC, pseq["seq"]).group(0))
            start = start + start_adj
            stop = stop + start_adj - stop_adj
            parsed[seq_type]["seq"] = pseq["seq"][start:stop]
    if len(parsed["query"]["seq"]) != len(parsed["hit"]["seq"]):
        raise ValueError(
            "Length mismatch: %r %r"
            % (len(parsed["query"]["seq"]), len(parsed["hit"]["seq"]))
        )
    if "similarity" in hsp.aln_annotation:
        
        
        hsp.aln_annotation["similarity"] = hsp.aln_annotation["similarity"][start:]
        
        assert len(hsp.aln_annotation["similarity"]) == len(parsed["hit"]["seq"])

    
    assert parsed["query"]["_type"] == parsed["hit"]["_type"]
    type_val = parsed["query"]["_type"]  
    molecule_type = "DNA" if type_val == "D" else "protein"
    setattr(hsp.fragment, "molecule_type", molecule_type)

    for seq_type in ("hit", "query"):
        
        start = int(parsed[seq_type]["_start"])
        end = int(parsed[seq_type]["_stop"])

        setattr(hsp.fragment, seq_type + "_start", min(start, end) - 1)
        setattr(hsp.fragment, seq_type + "_end", max(start, end))
        
        setattr(hsp.fragment, seq_type, parsed[seq_type]["seq"])

        if molecule_type != "protein":
            
            
            if start <= end:
                setattr(hsp.fragment, seq_type + "_strand", 1)
            else:
                setattr(hsp.fragment, seq_type + "_strand", -1)
        else:
            setattr(hsp.fragment, seq_type + "_strand", 0)


def _get_aln_slice_coords(parsed_hsp):
    """Get HSPs sequences (PRIVATE).

    To get the actual pairwise alignment sequences, we must first
    translate the un-gapped sequence based coordinates into positions
    in the gapped sequence (which may have a flanking region shown
    using leading - characters).  To date, I have never seen any
    trailing flanking region shown in the m10 file, but the
    following code should also cope with that.

    Note that this code seems to work fine even when the "sq_offset"
    entries are present as a result of using the -X command line option.
    """
    seq = parsed_hsp["seq"]
    seq_stripped = seq.strip("-")
    disp_start = int(parsed_hsp["_display_start"])
    start = int(parsed_hsp["_start"])
    stop = int(parsed_hsp["_stop"])

    if start <= stop:
        start = start - disp_start
        stop = stop - disp_start + 1
    else:
        start = disp_start - start
        stop = disp_start - stop + 1
    stop += seq_stripped.count("-")
    if not (0 <= start and start < stop and stop <= len(seq_stripped)):
        raise ValueError(
            "Problem with sequence start/stop,\n%s[%i:%i]\n%s"
            % (seq, start, stop, parsed_hsp)
        )
    return start, stop


class FastaM10Parser:
    """Parser for Bill Pearson's FASTA suite's -m 10 output."""

    def __init__(self, handle, __parse_hit_table=False):
        """Initialize the class."""
        self.handle = handle
        self._preamble = self._parse_preamble()

    def __iter__(self):
        """Iterate over FastaM10Parser object yields query results."""
        for qresult in self._parse_qresult():
            
            qresult.description = qresult.description
            yield qresult

    def _parse_preamble(self):
        """Parse the Fasta preamble for Fasta flavor and version (PRIVATE)."""
        preamble = {}
        while True:
            line = self.handle.readline()
            
            if line.startswith("Query"):
                break
            
            elif line.startswith(" version"):
                preamble["version"] = line.split(" ")[2]
            else:
                
                flav_match = re.match(_RE_FLAVS, line.lower())
                if flav_match:
                    preamble["program"] = flav_match.group(0)
        self.line = line

        return preamble

    def __parse_hit_table(self):
        """Parse hit table rows."""
        
        hit_rows = []
        while True:
            line = self.handle.readline()
            if (not line) or line.strip():
                break
            hit_rows.append("")
        self.line = line
        return hit_rows

    def _parse_qresult(self):
        """Parse query result (PRIVATE)."""
        
        qresult = None
        hit_rows = []
        
        state_QRES_NEW = 1
        state_QRES_HITTAB = 3
        state_QRES_CONTENT = 5
        state_QRES_END = 7

        line = self.line

        while True:
            
            if line.startswith("The best scores are:"):
                qres_state = state_QRES_HITTAB
            
            elif line.strip() == ">>>///" or not line:
                qres_state = state_QRES_END
            
            elif not line.startswith(">>>") and ">>>" in line:
                qres_state = state_QRES_NEW
            
            elif line.startswith(">>>") and line.strip() != ">>><<<":
                qres_state = state_QRES_CONTENT
            
            else:
                qres_state = None

            if qres_state is not None:
                if qres_state == state_QRES_HITTAB:
                    
                    hit_rows = self.__parse_hit_table()
                    line = self.handle.readline()

                elif qres_state == state_QRES_END:
                    yield _set_qresult_hits(qresult, hit_rows)
                    break

                elif qres_state == state_QRES_NEW:
                    
                    if qresult is not None:
                        yield _set_qresult_hits(qresult, hit_rows)
                    regx = re.search(_RE_ID_DESC_SEQLEN, line)
                    query_id = regx.group(1)
                    seq_len = regx.group(3)
                    desc = regx.group(2)
                    qresult = QueryResult(id=query_id)
                    qresult.seq_len = int(seq_len)
                    
                    line = self.handle.readline()
                    qresult.target = [x for x in line.split(" ") if x][1].strip()
                    if desc is not None:
                        qresult.description = desc
                    
                    for key, value in self._preamble.items():
                        setattr(qresult, key, value)
                    line = self.handle.readline()

                elif qres_state == state_QRES_CONTENT:
                    assert line[3:].startswith(qresult.id), line
                    for hit, strand in self._parse_hit(query_id):
                        
                        hit.description = hit.description
                        hit.query_description = qresult.description
                        
                        if hit.id not in qresult:
                            qresult.append(hit)
                        
                        else:
                            
                            
                            for hsp in hit.hsps:
                                assert strand != hsp.query_strand
                                qresult[hit.id].append(hsp)
                    line = self.line

            else:
                line = self.handle.readline()

        self.line = line

    def _parse_hit(self, query_id):
        """Parse hit on query identifier (PRIVATE)."""
        while True:
            line = self.handle.readline()
            if line.startswith(">>"):
                break

        state = _STATE_NONE
        strand = None
        hsp_list = []
        hsp = None
        parsed_hsp = None
        hit_desc = None
        seq_len = None
        while True:
            
            
            self.line = self.handle.readline()
            if self.line.strip() in [">>><<<", ">>>///"] or (
                not self.line.startswith(">>>") and ">>>" in self.line
            ):
                
                if state == _STATE_HIT_BLOCK:
                    parsed_hsp["hit"]["seq"] += line.strip()
                elif state == _STATE_CONS_BLOCK:
                    hsp.aln_annotation["similarity"] += line.strip("\r\n")
                
                _set_hsp_seqs(hsp, parsed_hsp, self._preamble["program"])
                hit = Hit(hsp_list)
                hit.description = hit_desc
                hit.seq_len = seq_len
                yield hit, strand
                hsp_list = []
                break
            
            elif line.startswith(">>"):
                
                if hsp_list:
                    _set_hsp_seqs(hsp, parsed_hsp, self._preamble["program"])
                    hit = Hit(hsp_list)
                    hit.description = hit_desc
                    hit.seq_len = seq_len
                    yield hit, strand
                    hsp_list = []
                
                try:
                    hit_id, hit_desc = line[2:].strip().split(" ", 1)
                except ValueError:
                    hit_id = line[2:].strip().split(" ", 1)[0]
                    hit_desc = ""
                
                frag = HSPFragment(hit_id, query_id)
                hsp = HSP([frag])
                hsp_list.append(hsp)
                
                state = _STATE_NONE
                parsed_hsp = {"query": {}, "hit": {}}
            
            elif line.startswith(">--"):
                
                _set_hsp_seqs(hsp, parsed_hsp, self._preamble["program"])
                
                frag = HSPFragment(hit_id, query_id)
                hsp = HSP([frag])
                hsp_list.append(hsp)
                
                state = _STATE_NONE
                parsed_hsp = {"query": {}, "hit": {}}
            
            elif line.startswith(">"):
                if state == _STATE_NONE:
                    
                    if not query_id.startswith(line[1:].split(" ")[0]):
                        raise ValueError(f"{query_id!r} vs {line!r}")
                    state = _STATE_QUERY_BLOCK
                    parsed_hsp["query"]["seq"] = ""
                elif state == _STATE_QUERY_BLOCK:
                    
                    assert hit_id.startswith(line[1:].split(" ")[0])
                    state = _STATE_HIT_BLOCK
                    parsed_hsp["hit"]["seq"] = ""
            
            elif line.startswith("; al_cons"):
                state = _STATE_CONS_BLOCK
                hsp.fragment.aln_annotation["similarity"] = ""
            elif line.startswith(";"):
                
                
                
                regx = re.search(_RE_ATTR, line.strip())
                name = regx.group(1)
                value = regx.group(2)

                
                if state == _STATE_NONE:
                    if name in _HSP_ATTR_MAP:
                        attr_name, caster = _HSP_ATTR_MAP[name]
                        if caster is not str:
                            value = caster(value)
                        if name in ["_ident", "_sim"]:
                            value *= 100
                        setattr(hsp, attr_name, value)
                
                elif state == _STATE_QUERY_BLOCK:
                    parsed_hsp["query"][name] = value
                elif state == _STATE_HIT_BLOCK:
                    if name == "_len":
                        seq_len = int(value)
                    else:
                        parsed_hsp["hit"][name] = value
                
                else:
                    raise ValueError("Unexpected line: %r" % line)
            
            else:
                assert ">" not in line
                
                if state == _STATE_HIT_BLOCK:
                    parsed_hsp["hit"]["seq"] += line.strip()
                elif state == _STATE_QUERY_BLOCK:
                    parsed_hsp["query"]["seq"] += line.strip()
                elif state == _STATE_CONS_BLOCK:
                    hsp.fragment.aln_annotation["similarity"] += line.strip("\r\n")
                
                else:
                    raise ValueError("Unexpected line: %r" % line)
            line = self.line


class FastaM10Indexer(SearchIndexer):
    """Indexer class for Bill Pearson's FASTA suite's -m 10 output."""

    _parser = FastaM10Parser

    def __init__(self, filename):
        """Initialize the class."""
        SearchIndexer.__init__(self, filename)

    def __iter__(self):
        """Iterate over FastaM10Indexer; yields query results' keys, start offsets, offset lengths."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()
        qresult_key = None
        query_mark = b">>>"

        line = handle.readline()
        while True:
            end_offset = handle.tell()

            if not line.startswith(query_mark) and query_mark in line:
                regx = re.search(_RE_ID_DESC_SEQLEN_IDX, line)
                qresult_key = regx.group(1).decode()
                start_offset = end_offset - len(line)
            
            if qresult_key is not None:
                if not line:
                    yield qresult_key, start_offset, end_offset - start_offset
                    break
                line = handle.readline()
                if not line.startswith(query_mark) and query_mark in line:
                    yield qresult_key, start_offset, end_offset - start_offset
                    start_offset = end_offset
            else:
                line = handle.readline()

    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string."""
        handle = self._handle
        qresult_raw = b""
        query_mark = b">>>"

        
        handle.seek(0)
        line = handle.readline()
        while True:
            qresult_raw += line
            line = handle.readline()
            if not line.startswith(query_mark) and query_mark in line:
                break

        
        handle.seek(offset)
        line = handle.readline()
        while True:
            
            if not line:
                break
            qresult_raw += line

            line = handle.readline()
            
            if not line.startswith(query_mark) and query_mark in line:
                break

        
        return qresult_raw + b">>><<<\n"



if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
