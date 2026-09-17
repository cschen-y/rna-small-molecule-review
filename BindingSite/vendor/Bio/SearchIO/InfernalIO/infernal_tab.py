






"""Bio.SearchIO parser for Infernal tabular output format."""

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult

from ._base import _BaseInfernalParser
from Bio.File import _open_for_random_access

__all__ = ("InfernalTabParser", "InfernalTabIndexer")


_TAB_FORMAT = {
    1: (
        "target_name",
        "target_acc",
        "query_name",
        "query_acc",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "evalue",
        "inc",
        "description",
    ),
    2: (
        "idx",
        "target_name",
        "target_acc",
        "query_name",
        "query_acc",
        "clan",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "evalue",
        "inc",
        "olp",
        "anyidx",
        "afrct1",
        "afrct2",
        "winidx",
        "wfrct1",
        "wfrct2",
        "mdl_len",
        "seq_len",
        "description",
    ),
    3: (
        "target_name",
        "target_acc",
        "query_name",
        "query_acc",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "evalue",
        "inc",
        "mdl_len",
        "seq_len",
        "description",
    ),
}

_TABULAR_FMT_1_HEADER_FIELDS = [
    "target name",
    "accession",
    "query name",
    "accession",
    "mdl",
    "mdl from",
    "mdl to",
    "seq from",
    "seq to",
    "strand",
    "trunc",
    "pass",
    "gc",
    "bias",
    "score",
    "E-value",
    "inc",
    "description of target",
]
_TABULAR_FMT_2_HEADER_FIELDS = [
    "idx",
    "target name",
    "accession",
    "query name",
    "accession",
    "clan name",
    "mdl",
    "mdl from",
    "mdl to",
    "seq from",
    "seq to",
    "strand",
    "trunc",
    "pass",
    "gc",
    "bias",
    "score",
    "E-value",
    "inc",
    "olp",
    "anyidx",
    "afrct1",
    "afrct2",
    "winidx",
    "wfrct1",
    "wfrct2",
    "mdl len",
    "seq len",
    "description of target",
]
_TABULAR_FMT_3_HEADER_FIELDS = [
    "target name",
    "accession",
    "query name",
    "accession",
    "mdl",
    "mdl from",
    "mdl to",
    "seq from",
    "seq to",
    "strand",
    "trunc",
    "pass",
    "gc",
    "bias",
    "score",
    "E-value",
    "inc",
    "mdl len",
    "seq len",
    "description of target",
]


_COLUMN_QRESULT = {
    "query_name": ("id", str),
    "query_acc": ("accession", str),
    "seq_len": ("seq_len", int),
    "clan": ("clan", str),
}
_COLUMN_HIT = {
    "target_name": ("id", str),
    "target_acc": ("accession", str),
    "query_name": ("query_id", str),
    "description": ("description", str),
    "mdl_len": ("seq_len", int),
}
_COLUMN_HSP = {
    "score": ("bitscore", float),
    "evalue": ("evalue", float),
    "bias": ("bias", float),
    "gc": ("gc", float),
    "mdl": ("model", str),
    "trunc": ("truncated", str),
    "pass": ("pipeline_pass", int),
    "inc": ("is_included", lambda x: True if x == "!" else False),
    "olp": ("olp", str),
    "anyidx": ("anyidx", lambda x: None if x == "-" else int(x)),
    "afrct1": ("afrct1", lambda x: None if x == "-" else float(x)),
    "afrct2": ("afrct2", lambda x: None if x == "-" else float(x)),
    "winidx": ("winidx", lambda x: None if (x == "-" or x == '"') else int(x)),
    "wfrct1": ("wfrct1", lambda x: None if (x == "-" or x == '"') else float(x)),
    "wfrct2": ("wfrct2", lambda x: None if (x == "-" or x == '"') else float(x)),
}
_COLUMN_FRAG = {
    "mdl_from": ("query_start", int),
    "mdl_to": ("query_end", int),
    "seq_from": ("hit_start", int),
    "seq_to": ("hit_end", int),
    "strand": ("hit_strand", str),
}


def _infer_tabular_format(handle, is_byte=False):
    """Infer tabular format from the tabular file header (PRIVATE)."""
    
    
    
    
    
    
    
    
    

    def process_line(handle, is_byte):
        """Process a line from the handle, decoding byte stream."""
        line = handle.readline().strip()
        return line.decode() if is_byte else line

    header_label = process_line(handle, is_byte)
    header_delim = process_line(handle, is_byte)

    if header_label == "" or header_delim == "":
        raise ValueError(
            "Unexpected empty line in the header of Infernal tabular output."
        )

    
    
    indices = [idx for idx, char in enumerate(header_delim) if char == " "] + [
        len(header_delim)
    ]
    prev_idx = 1
    fields = []
    for idx in indices:
        fields.append(header_label[prev_idx:idx].strip())
        prev_idx = idx

    
    if fields == _TABULAR_FMT_1_HEADER_FIELDS:
        return 1
    elif fields == _TABULAR_FMT_2_HEADER_FIELDS:
        return 2
    elif fields == _TABULAR_FMT_3_HEADER_FIELDS:
        return 3
    else:
        raise ValueError(
            "Cannot determine the tabular format from the header. \
            The tabular file is likely incorrect or its format is unsupported. \
            Tabular format 1, 2 and 3 are currently supported)."
        )


class InfernalTabParser(_BaseInfernalParser):
    """Parser for the Infernal tabular format."""

    def __init__(self, handle, _fmt=None):
        """Initialize the class."""
        self.handle = handle
        if isinstance(_fmt, int):
            assert _fmt in {1, 2, 3}
            self.fmt = _fmt
        else:
            self.fmt = _infer_tabular_format(handle=handle)

    def __iter__(self):
        """Iterate over InfernalTabParser, yields query results."""
        
        self.line = self.handle.readline()
        while self.line.startswith("#"):
            self.line = self.handle.readline()
        
        if self.line:
            yield from self._parse_qresult()

    def _parse_row(self):
        """Return a dictionary of parsed row values (PRIVATE)."""
        cols = [x for x in self.line.strip().split(" ") if x]
        if len(cols) < len(_TAB_FORMAT[self.fmt]):
            raise ValueError(
                f"Less columns than expected for format {self.fmt}, only {len(cols)}"
            )
        
        cols[len(_TAB_FORMAT[self.fmt]) - 1] = " ".join(
            cols[len(_TAB_FORMAT[self.fmt]) - 1 :]
        )

        qresult, hit, hsp, frag = {}, {}, {}, {}
        for sname, value in zip(
            _TAB_FORMAT[self.fmt], cols[: len(_TAB_FORMAT[self.fmt])]
        ):
            
            
            for parsed_dict, mapping in (
                (qresult, _COLUMN_QRESULT),
                (hit, _COLUMN_HIT),
                (hsp, _COLUMN_HSP),
                (frag, _COLUMN_FRAG),
            ):
                
                if sname in mapping:
                    attr_name, caster = mapping[sname]
                    if caster is not str:
                        value = caster(value)
                    parsed_dict[attr_name] = value

        
        self._adjust_coords(frag)

        return {"qresult": qresult, "hit": hit, "hsp": hsp, "frag": frag}

    def _adjust_coords(self, frag):
        """Adjust start and end coordinates according to strand (PRIVATE)."""
        strand = frag["hit_strand"]
        assert strand is not None
        
        if strand == "-":
            hit_start = frag["hit_start"]
            hit_end = frag["hit_end"]
            frag["hit_start"] = hit_end
            frag["hit_end"] = hit_start
            frag["hit_strand"] = -1
        else:
            frag["hit_strand"] = 0

    def _parse_qresult(self):
        """Yield QueryResult objects (PRIVATE)."""
        
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        
        qres_state = None
        file_state = None
        cur_qid = None
        cur_hid = None
        
        prev_qid = None
        prev_hid = None
        
        cur, prev = None, None
        hit_dict = {}

        while True:
            
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            
            if self.line and not self.line.startswith("#"):
                cur = self._parse_row()
                cur_qid = cur["qresult"]["id"]
                cur_hid = cur["hit"]["id"]
            else:
                file_state = state_EOF
                
                cur_qid, cur_hid = None, None

            
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME

            
            
            if prev is not None:
                
                
                frag = HSPFragment(prev_hid, prev_qid)
                for attr, value in prev["frag"].items():
                    setattr(frag, attr, value)
                hsp = HSP([frag])
                for attr, value in prev["hsp"].items():
                    setattr(hsp, attr, value)
                self._add_hit_to_dict(prev["hit"], hsp, hit_dict)

                
                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    
                    qresult = QueryResult(self._hit_to_list(hit_dict), prev_qid)
                    for attr, value in prev["qresult"].items():
                        setattr(qresult, attr, value)
                    yield qresult
                    
                    if file_state == state_EOF:
                        break
                    hit_dict = {}

            self.line = self.handle.readline()


class InfernalTabIndexer(SearchIndexer):
    """Indexer class for Infernal tabular output."""

    _parser = InfernalTabParser

    def __init__(self, filename):
        """Initialize the class."""
        
        
        with _open_for_random_access(filename) as handle:
            fmt = _infer_tabular_format(handle, is_byte=True)
        self._query_id_idx = 3 if fmt == 2 else 2

        SearchIndexer.__init__(self, filename, _fmt=fmt)

    def __iter__(self):
        """Iterate over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        qresult_key = None
        start_offset = handle.tell()

        
        line = b"#"

        
        while line.startswith(b"#"):
            start_offset = handle.tell()
            line = handle.readline()

        
        while True:
            
            
            end_offset = handle.tell()

            if not line:
                break

            cols = [x for x in line.strip().split(b" ") if x]

            if qresult_key is None:
                qresult_key = cols[self._query_id_idx]
            else:
                curr_key = cols[self._query_id_idx]
                if curr_key != qresult_key:
                    adj_end = end_offset - len(line)
                    yield (qresult_key.decode(), start_offset, adj_end - start_offset)
                    qresult_key = curr_key
                    start_offset = adj_end

            line = handle.readline()
            if not line or line.startswith(b"#"):
                yield (qresult_key.decode(), start_offset, end_offset - start_offset)
                break

    def get_raw(self, offset):
        """Return the raw bytes string of a QueryResult object from the given offset."""
        handle = self._handle
        qresult_key = None
        qresult_raw = b""

        
        handle.seek(offset)
        while True:
            line = handle.readline()
            if not line or line.startswith(b"#"):
                break
            cols = [x for x in line.strip().split(b" ") if x]
            if qresult_key is None:
                qresult_key = cols[self._query_id_idx]
            else:
                curr_key = cols[self._query_id_idx]
                if curr_key != qresult_key:
                    break
            qresult_raw += line

        return qresult_raw



if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
