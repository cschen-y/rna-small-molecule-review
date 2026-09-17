







"""Bio.SearchIO parser for Infernal plain text output format."""

import operator
import re

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult
from Bio.SearchIO._utils import read_forward

from ._base import _BaseInfernalParser

__all__ = ("InfernalTextParser", "InfernalTextIndexer")



_RE_PROGRAM = re.compile(r"^# .*?(\w?cm\w+) :: .*$")

_RE_VERSION = re.compile(r"# INFERNAL+ ([\w+\.]{2,}) .*$")

_RE_OPT = re.compile(r"^# (.+):\s+(.+)$")

_RE_NUMERIC = re.compile(r"\d+")
_RE_NOT_NUMERIC = re.compile(r"\D")

_RE_LETTERS = re.compile(r"[^A-Za-z]")

_RE_SPLIT_ALN = re.compile(r"(\*\[[ 0-9]+\]\*)+")





_DIV_HEADER_OPT = "# - - -"

_DIV_QUERY_END = "//"
_DIV_QUERY_START = "Query:"

_DIV_HITS_END_CM = "Internal CM pipeline statistics summary:"
_DIV_HITS_END_HMM = "Internal HMM-only pipeline statistics summary:"
_DIV_HIT_SCORE_TABLE = "Hit scores:"
_DIV_HIT_ALIGNMENT = "Hit alignments:"
_DIV_NO_HIT = "   [No hits detected that satisfy reporting thresholds]"
_DIV_TABLE_START = " ----   --------- ------"
_DIV_INC_THRESHOLD = " ------ inclusion threshold ------"

_DIV_ALIGNMENT_START = ">> "


class InfernalTextParser(_BaseInfernalParser):
    """Parser for the Infernal text output."""

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        self.line = read_forward(self.handle)
        self._meta = self._parse_header()

    def __iter__(self):
        """Iterate over query results."""
        yield from self._parse_qresult()

    def _read_until(self, bool_func):
        """Read the file handle until the given function returns True (PRIVATE)."""
        while True:
            if not self.line or bool_func(self.line):
                return
            self.line = read_forward(self.handle)

    def _parse_header(self):
        """Parse Infernal header (PRIVATE)."""
        meta = {}
        
        
        meta["show alignments in output"] = "yes"
        
        
        has_opts = False
        while True:
            
            if not self.line.startswith("#"):
                break
            
            
            elif self.line.startswith(_DIV_HEADER_OPT):
                if not has_opts:
                    
                    
                    has_opts = True
                else:
                    
                    
                    break
            elif not has_opts:
                
                regx = re.search(_RE_PROGRAM, self.line)
                if regx:
                    meta["program"] = regx.group(1)
                
                regx = re.search(_RE_VERSION, self.line)
                if regx:
                    meta["version"] = regx.group(1)
            elif has_opts:
                regx = re.search(_RE_OPT, self.line)
                
                if "target" in regx.group(1):
                    meta["target"] = regx.group(2).strip()
                else:
                    meta[regx.group(1)] = regx.group(2)

            self.line = read_forward(self.handle)

        return meta

    def _parse_qresult(self):
        """Parse a Infernal query block (PRIVATE)."""
        self._read_until(lambda line: line.startswith(_DIV_QUERY_START))

        while self.line:
            
            if self.line.startswith(_DIV_QUERY_START):
                qid = self.line.strip().split()[1]
                qlen = int(re.sub(_RE_NOT_NUMERIC, "", self.line.strip().split()[-1]))

                
                qresult_attrs = {
                    "id": qid,
                    "seq_len": qlen,
                    "program": self._meta.get("program"),
                    "version": self._meta.get("version"),
                    "target": self._meta.get("target"),
                }
            else:
                self.line = read_forward(self.handle)

            
            qdesc = "<unknown description>"  
            while not self.line.startswith(_DIV_HIT_SCORE_TABLE):
                self.line = read_forward(self.handle)

                if self.line.startswith("Accession:"):
                    acc = self.line.strip().split(" ", 1)[1]
                    qresult_attrs["accession"] = acc.strip()
                elif self.line.startswith("Description:"):
                    qdesc = self.line.strip().split(" ", 1)[1].strip()
                    qresult_attrs["description"] = qdesc

            
            
            hit_list = []
            while self.line and not self.line.startswith(_DIV_QUERY_END):
                hit_list = self._parse_hit(qid, qdesc)
                
                if self.line.startswith(_DIV_HITS_END_CM) or self.line.startswith(
                    _DIV_HITS_END_HMM
                ):
                    while self.line and not self.line.startswith(_DIV_QUERY_END):
                        self.line = read_forward(self.handle)
            
            qresult = QueryResult(id=qid, hits=hit_list)
            for attr, value in qresult_attrs.items():
                setattr(qresult, attr, value)
            yield qresult
            self.line = read_forward(self.handle)

            
            
            if "[ok]" in self.line:
                break

    def _parse_hit(self, qid, qdesc):
        """Parse an Infernal hit (PRIVATE)."""
        
        hit_end = False
        in_score_table = False
        parsing_hits = False
        prev_line = None
        
        
        
        
        hit_dict = {}

        
        if self._meta["show alignments in output"] == "no":
            div_hit_start = _DIV_HIT_SCORE_TABLE
        else:
            div_hit_start = _DIV_ALIGNMENT_START

        while True:
            if not self.line:
                raise ValueError("Unexpected end of file")
            
            elif self.line.startswith(_DIV_NO_HIT):
                while True:
                    self.line = read_forward(self.handle)
                    if self.line.startswith(_DIV_HITS_END_CM) or self.line.startswith(
                        _DIV_HITS_END_HMM
                    ):
                        hit_end = True
                        return []
            
            elif self.line.startswith(div_hit_start):
                
                if self._meta["show alignments in output"] == "no":
                    assert not in_score_table
                    self._read_until(lambda line: line.startswith(_DIV_TABLE_START))
                    self.line = read_forward(self.handle)
                    parsing_hits = in_score_table = True
                else:
                    self._parse_hit_from_alignment(qid, hit_dict)
                    parsing_hits = True
            
            elif self.line.startswith(_DIV_HITS_END_CM) or self.line.startswith(
                _DIV_HITS_END_HMM
            ):
                hit_end = True
                parsing_hits = False
            
            
            
            else:
                
                if in_score_table:
                    if self.line.strip():
                        prev_line = self.line
                    else:
                        parsing_hits = in_score_table = False

                
                self.line = self.handle.readline()

            
            if in_score_table and prev_line is not None:
                if not prev_line.startswith(_DIV_INC_THRESHOLD):
                    self._parse_scores_table_row(prev_line, qid, hit_dict)

            
            
            if hit_end:
                return self._hit_to_list(hit_dict)

    def _parse_hit_from_alignment(self, qid, hit_dict):
        """Parse an Infernal hit alignment (PRIVATE)."""
        hid, hdesc = self.line[len(_DIV_ALIGNMENT_START) :].split("  ", 1)
        hdesc = hdesc.strip()

        
        self._read_until(lambda line: line.startswith(_DIV_TABLE_START))
        self.line = read_forward(self.handle)

        
        row = [x for x in self.line.strip().split() if x]
        assert len(row) == 16

        
        hit_attrs = {"id": hid, "query_id": qid, "description": hdesc}
        hsp_attrs = {
            "evalue": float(row[2]),
            "bitscore": float(row[3]),
            "bias": float(row[4]),
            "model": row[5],
            "truncated": row[14],
            "gc": float(row[15]),
            "avg_acc": float(row[13]),
            "query_endtype": row[8],
            "hit_endtype": row[12],
            "is_included": True if row[1] == "!" else False,
        }
        query_start = int(row[6])
        query_end = int(row[7])
        hit_start = int(row[9]) if row[11] == "+" else int(row[10])
        hit_end = int(row[10]) if row[11] == "+" else int(row[9])
        hit_strand = 0 if row[11] == "+" else -1

        
        self.line = read_forward(self.handle)

        
        frag_list = self._parse_aln_block(
            hit_attrs["id"],
            hit_attrs["query_id"],
            hsp_attrs["model"],
            query_start,
            query_end,
            hit_start,
            hit_end,
            hit_strand,
        )
        hsp = HSP(frag_list)
        for attr, value in hsp_attrs.items():
            setattr(hsp, attr, value)

        
        self._add_hit_to_dict(hit_attrs, hsp, hit_dict)

    def _parse_scores_table_row(self, row, qid, hit_dict):
        """Parse an Infernal hit scores table (when used with --noali) (PRIVATE)."""

        
        row = [x for x in row.strip().split(" ") if x]
        
        if len(row) > 12:
            row[12] = " ".join(row[12:])
        
        elif len(row) < 12:
            row.append("")
            assert len(row) == 12

        
        hit_attrs = {"id": row[5], "query_id": qid, "description": row[12]}
        hsp_attrs = {
            "evalue": float(row[2]),
            "bitscore": float(row[3]),
            "bias": float(row[4]),
            "model": row[9],
            "truncated": row[10],
            "gc": float(row[11]),
            "is_included": True if row[1] == "!" else False,
        }
        hsp_frag_attrs = {
            "hit_start": int(row[6]) if row[8] == "+" else int(row[7]),
            "hit_end": int(row[7]) if row[8] == "+" else int(row[6]),
            "hit_strand": 0 if row[8] == "+" else -1,
        }

        
        hsp_frag = HSPFragment(row[5], qid)
        for attr, value in hsp_frag_attrs.items():
            setattr(hsp_frag, attr, value)

        
        hsp = HSP([hsp_frag])
        for attr, value in hsp_attrs.items():
            setattr(hsp, attr, value)

        
        self._add_hit_to_dict(hit_attrs, hsp, hit_dict)

    def _parse_aln_block(
        self, hid, qid, model, query_start, query_end, hit_start, hit_end, hit_strand
    ):
        """Parse a Infernal HSP alignment block (PRIVATE)."""
        frag_list = []
        model_seq = ""
        hit_seq = ""
        if model == "cm":
            annot = {"NC": "", "CS": "", "similarity": "", "PP": ""}
        else:
            annot = {"CS": "", "similarity": "", "PP": ""}
        while True:
            
            if (
                self.line.startswith(_DIV_ALIGNMENT_START)
                or self.line.startswith(_DIV_HITS_END_CM)
                or self.line.startswith(_DIV_HITS_END_HMM)
            ):
                
                
                
                
                

                
                
                local_aln_idx = [
                    (0, 0)
                ]  
                local_aln_idx += [
                    (m.start(0), m.end(0))
                    for m in re.finditer(_RE_SPLIT_ALN, model_seq)
                ]

                prev_hit_start = hit_start if hit_strand == 0 else hit_end
                prev_model_start = query_start
                hsps = []

                for i in range(len(local_aln_idx)):
                    local_start = local_aln_idx[i][1]
                    local_end = (
                        local_aln_idx[i + 1][0] if i + 1 < len(local_aln_idx) else None
                    )

                    
                    op = operator.add if hit_strand == 0 else operator.sub
                    cur_hit_seq = hit_seq[local_start:local_end]
                    cur_hit_gap_size = self._local_aln_gap_size(
                        local_aln_idx[i], hit_seq
                    )
                    cur_hit_start = op(prev_hit_start, cur_hit_gap_size)
                    cur_hit_end = op(
                        cur_hit_start, len(re.sub(_RE_LETTERS, "", cur_hit_seq))
                    )
                    
                    if hit_strand == 0 and i == 0:
                        cur_hit_end -= 1
                    if hit_strand == -1 and i == 0:
                        cur_hit_end += 1
                    prev_hit_start = cur_hit_end
                    
                    cur_model_seq = model_seq[local_start:local_end].replace(".", "-")
                    cur_model_gap_size = self._local_aln_gap_size(
                        local_aln_idx[i], model_seq
                    )
                    cur_model_start = prev_model_start + cur_model_gap_size
                    cur_model_end = cur_model_start + len(
                        re.sub(_RE_LETTERS, "", cur_model_seq)
                    )
                    if i == 0:
                        cur_model_end -= 1
                    prev_model_start = cur_model_end
                    
                    cur_annot = {k: v[local_start:local_end] for k, v in annot.items()}

                    
                    frag = HSPFragment(hid, qid)
                    frag.query = cur_model_seq
                    frag.hit = cur_hit_seq
                    frag.query_start = cur_model_start
                    frag.query_end = cur_model_end
                    frag.hit_start = cur_hit_start if hit_strand == 0 else cur_hit_end
                    frag.hit_end = cur_hit_end if hit_strand == 0 else cur_hit_start
                    frag.hit_strand = hit_strand
                    frag.aln_annotation = cur_annot

                    frag_list.append(frag)

                return frag_list

            
            
            block_size = 6 if model == "cm" else 5
            offset = 1 if model == "cm" else 0  
            lines = [None] * block_size
            for i in range(block_size):
                lines[i] = self.line
                self.line = read_forward(self.handle)

            
            blklen = len(lines[4 + offset].strip().split()[0])
            blkstart = len(lines[4 + offset]) - blklen - 4
            blkend = len(lines[4 + offset]) - 4
            model_seq += lines[1 + offset][blkstart:blkend]
            hit_seq += lines[3 + offset][blkstart:blkend]
            
            if model == "cm":
                annot["NC"] += lines[0][blkstart:blkend]
            annot["CS"] += lines[0 + offset][blkstart:blkend]
            annot["similarity"] += lines[2 + offset][blkstart:blkend]
            annot["PP"] += lines[4 + offset][blkstart:blkend]

    def _local_aln_gap_size(self, cur_aln_idx, seq):
        """Calculate the gap size between the local alignments (PRIVATE)."""
        gap_len = 0
        if cur_aln_idx[1] > 0:
            gap_len = sum(
                [
                    int(n)
                    for n in re.findall(
                        _RE_NUMERIC, seq[cur_aln_idx[0] : cur_aln_idx[1]]
                    )
                ]
            )
            assert gap_len > 0
        return gap_len


class InfernalTextIndexer(SearchIndexer):
    """Indexer class for Infernal plain text output."""

    _parser = InfernalTextParser

    def __init__(self, *args, **kwargs):
        """Initialize the class."""
        super().__init__(*args, **kwargs)
        self._preamble = b""

    def __iter__(self):
        """Iterate over InfernalTextIndexer; yields query results' key, offsets, 0."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()

        while True:
            line = read_forward(handle)
            end_offset = handle.tell()

            if line.startswith(_DIV_QUERY_START.encode()):
                qresult_key = line.strip().split()[1]
                
                
                start_offset = end_offset - len(line)
            elif line.startswith(_DIV_QUERY_END.encode()):
                yield qresult_key.decode(), start_offset, 0
                start_offset = end_offset
            elif not line:
                break

    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string."""
        handle = self._handle
        qresult_raw = b""

        
        if not self._preamble:
            handle.seek(0)
            while True:
                line = handle.readline()
                if line.startswith(_DIV_QUERY_START.encode()):
                    break
                self._preamble += line

        qresult_raw += self._preamble

        
        handle.seek(offset)
        while True:
            
            line = handle.readline()
            qresult_raw += line

            
            if line.startswith(_DIV_QUERY_END.encode()) or not line:
                break

        return qresult_raw



if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
