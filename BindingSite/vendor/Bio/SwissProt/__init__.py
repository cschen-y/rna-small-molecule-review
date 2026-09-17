



"""Code to work with the sprotXX.dat file from SwissProt.

https://web.expasy.org/docs/userman.html

Classes:
 - Record             Holds SwissProt data.
 - Reference          Holds reference data from a SwissProt record.

Functions:
 - read               Read one SwissProt record
 - parse              Read multiple SwissProt records

"""

import io
import re

from Bio.SeqFeature import Position
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation


class SwissProtParserError(ValueError):
    """An error occurred while parsing a SwissProt file."""

    def __init__(self, *args, line=None):
        """Create a SwissProtParserError object with the offending line."""
        super().__init__(*args)
        self.line = line


class Record:
    """Holds information from a SwissProt record.

    Attributes:
     - entry_name        Name of this entry, e.g. RL1_ECOLI.
     - data_class        Either 'STANDARD' or 'PRELIMINARY'.
     - molecule_type     Type of molecule, 'PRT',
     - sequence_length   Number of residues.
     - accessions        List of the accession numbers, e.g. ['P00321']
     - created           A tuple of (date, release).
     - sequence_update   A tuple of (date, release).
     - annotation_update A tuple of (date, release).
     - description       Free-format description.
     - gene_name         A list of dictionaries with keys 'Name', 'Synonyms',
                         'OrderedLocusNames' and 'ORFNames'.
     - organism          The source of the sequence.
     - organelle         The origin of the sequence.
     - organism_classification  The taxonomy classification.  List of strings.
       (http://www.ncbi.nlm.nih.gov/Taxonomy/)
     - taxonomy_id       A list of NCBI taxonomy id's.
     - host_organism     A list of names of the hosts of a virus, if any.
     - host_taxonomy_id  A list of NCBI taxonomy id's of the hosts, if any.
     - references        List of Reference objects.
     - comments          List of strings.
     - cross_references  List of tuples (db, id1[, id2][, id3]).  See the docs.
     - keywords          List of the keywords.
     - features          List of tuples (key name, from, to, description).
       from and to can be either integers for the residue
       numbers, '<', '>', or '?'
     - protein_existence Numerical value describing the evidence for the existence of the protein.
     - seqinfo           tuple of (length, molecular weight, CRC32 value)
     - sequence          The sequence.

    Examples
    --------
    >>> from Bio import SwissProt
    >>> example_filename = "SwissProt/P68308.txt"
    >>> with open(example_filename) as handle:
    ...     records = SwissProt.parse(handle)
    ...     for record in records:
    ...         print(record.entry_name)
    ...         print(record.accessions)
    ...         print(record.keywords)
    ...         print(record.organism)
    ...         print(record.sequence[:20] + "...")
    ...
    NU3M_BALPH
    ['P68308', 'P24973']
    ['Electron transport', 'Membrane', 'Mitochondrion', 'Mitochondrion inner membrane', 'NAD', 'Respiratory chain', 'Translocase', 'Transmembrane', 'Transmembrane helix', 'Transport', 'Ubiquinone']
    Balaenoptera physalus (Fin whale) (Balaena physalus).
    MNLLLTLLTNTTLALLLVFI...

    """

    def __init__(self):
        """Initialize the class."""
        self.entry_name = None
        self.data_class = None
        self.molecule_type = None
        self.sequence_length = None

        self.accessions = []
        self.created = None
        self.sequence_update = None
        self.annotation_update = None

        self.description = []
        self.gene_name = []
        self.organism = []
        self.organelle = ""
        self.organism_classification = []
        self.taxonomy_id = []
        self.host_organism = []
        self.host_taxonomy_id = []
        self.references = []
        self.comments = []
        self.cross_references = []
        self.keywords = []
        self.features = []
        self.protein_existence = ""

        self.seqinfo = None
        self.sequence = ""


class Reference:
    """Holds information from one reference in a SwissProt entry.

    Attributes:
     - number      Number of reference in an entry.
     - evidence    Evidence code.  List of strings.
     - positions   Describes extent of work.  List of strings.
     - comments    Comments.  List of (token, text).
     - references  References.  List of (dbname, identifier).
     - authors     The authors of the work.
     - title       Title of the work.
     - location    A citation for the work.

    """

    def __init__(self):
        """Initialize the class."""
        self.number = None
        self.positions = []
        self.comments = []
        self.references = []
        self.authors = []
        self.title = []
        self.location = []


class FeatureTable(SeqFeature):
    """Stores feature annotations for specific regions of the sequence.

    This is a subclass of SeqFeature, defined in Bio.SeqFeature, where the
    attributes are used as follows:

     - ``location``: location of the feature on the canonical or isoform
       sequence; the location is stored as an instance of SimpleLocation,
       defined in Bio.SeqFeature, with the ref attribute set to the isoform
       ID referring to the canonical or isoform sequence on which the feature
       is defined
     - ``id``: unique and stable identifier (FTId), only provided for features
       belonging to the types CARBOHYD, CHAIN, PEPTIDE, PROPEP, VARIANT, or
       VAR_SEQ
     - ``type``: indicates the type of feature, as defined by the UniProt
       Knowledgebase documentation:

        - ACT_SITE: amino acid(s) involved in the activity of an enzyme
        - BINDING:  binding site for any chemical group
        - CARBOHYD: glycosylation site; an FTId identifier to the GlyConnect
          database is provided if annotated there
        - CA_BIND:  calcium-binding region
        - CHAIN:    polypeptide chain in the mature protein
        - COILED:   coiled-coil region
        - COMPBIAS: compositionally biased region
        - CONFLICT: different sources report differing sequences
        - CROSSLNK: posttransationally formed amino acid bond
        - DISULFID: disulfide bond
        - DNA_BIND: DNA-binding region
        - DOMAIN:   domain, defined as a specific combination of secondary
          structures organized into a characteristic three-dimensional
          structure or fold
        - INIT_MET: initiator methionine
        - INTRAMEM: region located in a membrane without crossing it
        - HELIX:    alpha-, 3(10)-, or pi-helix secondary structure
        - LIPID:    covalent binding of a lipid moiety
        - METAL:    binding site for a metal ion
        - MOD_RES:  posttranslational modification (PTM) of a residue,
          annotated by the controlled vocabulary defined by the ptmlist.txt
          document on the UniProt website
        - MOTIF:    short sequence motif of biological interest
        - MUTAGEN:  site experimentally altered by mutagenesis
        - NON_CONS: non-consecutive residues
        - NON_STD:  non-standard amino acid
        - NON_TER:  the residue at an extremity of the sequence is not the
          terminal residue
        - NP_BIND:  nucleotide phosphate-binding region
        - PEPTIDE:  released active mature polypeptide
        - PROPEP:   any processed propeptide
        - REGION:   region of interest in the sequence
        - REPEAT:   internal sequence repetition
        - SIGNAL:   signal sequence (prepeptide)
        - SITE:     amino-acid site of interest not represented by another
          feature key
        - STRAND:   beta-strand secondary structure; either a hydrogen-bonded
          extended beta strand or a residue in an isolated beta-bridge
        - TOPO_DOM: topological domain
        - TRANSIT:  transit peptide (mitochondrion, chloroplast, thylakoid,
          cyanelle, peroxisome, etc.)
        - TRANSMEM: transmembrane region
        - TURN:     H-bonded turn (3-, 4-, or 5-turn)
        - UNSURE:   uncertainties in the sequence
        - VARIANT:  sequence variant; an FTId is provided for protein sequence
          variants of Hominidae (great apes and humans)
        - VAR_SEQ:  sequence variant produced by alternative splicing,
          alternative promoter usage, alternative initiation, or ribosomal
          frameshifting
        - ZN_FING:  zinc finger region

     - ``qualifiers``: a dictionary of additional information, which may include
       the feature evidence and free-text notes. While SwissProt includes the
       feature identifier code (FTId) as a qualifier, it is stored as the
       attribute ID of the FeatureTable object.

    """


def parse(source):
    """Read multiple SwissProt records from file.

    Argument source is a file-like object or a path to a file.

    Returns a generator object which yields Bio.SwissProt.Record() objects.
    """
    handle = _open(source)
    try:
        while True:
            record = _read(handle)
            if not record:
                return
            yield record
    finally:
        if handle is not source:
            handle.close()


def read(source):
    """Read one SwissProt record from file.

    Argument source is a file-like object or a path to a file.

    Returns a Record() object.
    """
    handle = _open(source)
    try:
        record = _read(handle)
        if not record:
            raise ValueError("No SwissProt record found")
        
        
        try:
            next(handle)
        except StopIteration:
            return record
        raise ValueError("More than one SwissProt record found")
    finally:
        if handle is not source:
            handle.close()





def _open(source):
    try:
        handle = open(source)
        return handle
    except TypeError:
        handle = source
        if handle.read(0) == "":
            
            return handle
        
        return io.TextIOWrapper(handle, encoding="ASCII")


def _read(handle):
    record = None
    unread = ""
    try:
        line = next(handle)
    except StopIteration:
        return record
    key, value = line[:2], line[5:].rstrip()
    if key != "ID":
        raise SwissProtParserError("Failed to find ID in first line", line=line)
    record = Record()
    _read_id(record, line)
    _sequence_lines = []
    for line in handle:
        key, value = line[:2], line[5:].rstrip()
        if unread:
            value = unread + " " + value
            unread = ""
        if key == "AC":
            accessions = value.rstrip(";").split("; ")
            record.accessions.extend(accessions)
        elif key == "DT":
            _read_dt(record, line)
        elif key == "DE":
            record.description.append(value.strip())
        elif key == "GN":
            if value == "and":
                record.gene_name.append("")
            else:
                if len(record.gene_name) == 0:
                    record.gene_name.append("")
                record.gene_name[-1] += value + " "
        elif key == "OS":
            record.organism.append(value)
        elif key == "OG":
            record.organelle += line[5:]
        elif key == "OC":
            cols = value.rstrip(";.").split("; ")
            record.organism_classification.extend(cols)
        elif key == "OX":
            _read_ox(record, line)
        elif key == "OH":
            _read_oh(record, line)
        elif key == "RN":
            reference = Reference()
            _read_rn(reference, value)
            record.references.append(reference)
        elif key == "RP":
            assert record.references, "RP: missing RN"
            record.references[-1].positions.append(value)
        elif key == "RC":
            assert record.references, "RC: missing RN"
            reference = record.references[-1]
            unread = _read_rc(reference, value)
        elif key == "RX":
            assert record.references, "RX: missing RN"
            reference = record.references[-1]
            _read_rx(reference, value)
        elif key == "RL":
            assert record.references, "RL: missing RN"
            reference = record.references[-1]
            reference.location.append(value)
        
        
        
        elif key == "RA":
            assert record.references, "RA: missing RN"
            reference = record.references[-1]
            reference.authors.append(value)
        elif key == "RG":
            assert record.references, "RG: missing RN"
            reference = record.references[-1]
            reference.authors.append(value)
        elif key == "RT":
            assert record.references, "RT: missing RN"
            reference = record.references[-1]
            reference.title.append(value)
        elif key == "CC":
            _read_cc(record, line)
        elif key == "DR":
            _read_dr(record, value)
        elif key == "PE":
            _read_pe(record, value)
        elif key == "KW":
            _read_kw(record, value)
        elif key == "FT":
            _read_ft(record, line)
        elif key == "SQ":
            cols = value.split()
            assert len(cols) == 7, f"I don't understand SQ line {line}"
            
            record.seqinfo = int(cols[1]), int(cols[3]), cols[5]
        elif key == "  ":
            _sequence_lines.append(value.replace(" ", "").rstrip())
        elif key == "//":
            _read_gn(record)
            
            record.description = " ".join(record.description)
            record.organism = " ".join(record.organism)
            record.organelle = record.organelle.rstrip()
            for reference in record.references:
                reference.authors = " ".join(reference.authors).rstrip(";")
                if reference.title:
                    title = reference.title[0]
                    for fragment in reference.title[1:]:
                        if not title.endswith("-"):
                            title += " "
                        title += fragment
                    title = title.rstrip(";")
                    if title.startswith('"') and title.endswith('"'):
                        title = title[1:-1]  
                else:
                    title = ""
                reference.title = title
                reference.location = " ".join(reference.location)
            record.sequence = "".join(_sequence_lines)
            return record
        elif key == "**":
            
            
            
            
            
            
            
            pass
        else:
            raise SwissProtParserError(f"Unknown keyword {key!r} found", line=line)
    if record:
        raise ValueError("Unexpected end of stream.")


def _read_gn(record):
    for i, text in enumerate(record.gene_name):
        tokens = text.rstrip("; ").split("; ")
        gene_name = {}
        for token in tokens:
            
            
            key, value = token.strip().split("=", 1)
            if key == "Name":
                gene_name["Name"] = value
            else:
                assert key in ("Synonyms", "OrderedLocusNames", "ORFNames")
                gene_name[key] = value.split(", ")
        record.gene_name[i] = gene_name


def _read_id(record, line):
    cols = line[5:].split()
    
    
    
    
    
    if len(cols) == 5:
        record.entry_name = cols[0]
        record.data_class = cols[1].rstrip(";")
        record.molecule_type = cols[2].rstrip(";")
        record.sequence_length = int(cols[3])
    elif len(cols) == 4:
        record.entry_name = cols[0]
        record.data_class = cols[1].rstrip(";")
        record.molecule_type = None
        record.sequence_length = int(cols[2])
    else:
        raise SwissProtParserError("ID line has unrecognised format", line=line)
    
    allowed = ("STANDARD", "PRELIMINARY", "IPI", "Reviewed", "Unreviewed")
    if record.data_class not in allowed:
        message = f"Unrecognized data class {record.data_class!r}"
        raise SwissProtParserError(message, line=line)

    
    
    if record.molecule_type not in (None, "PRT"):
        message = f"Unrecognized molecule type {record.molecule_type!r}"
        raise SwissProtParserError(message, line=line)


def _read_dt(record, line):
    value = line[5:]
    uprline = value.upper()
    cols = value.rstrip().split()
    if (
        "CREATED" in uprline
        or "LAST SEQUENCE UPDATE" in uprline
        or "LAST ANNOTATION UPDATE" in uprline
    ):
        
        
        
        
        
        
        
        
        
        

        
        
        
        
        uprcols = uprline.split()
        rel_index = -1
        for index in range(len(uprcols)):
            if "REL." in uprcols[index]:
                rel_index = index
        assert rel_index >= 0, f"Could not find Rel. in DT line: {line}"
        version_index = rel_index + 1
        
        str_version = cols[version_index].rstrip(",")
        
        if str_version == "":
            version = 0
        
        elif "." in str_version:
            version = str_version
        
        else:
            version = int(str_version)
        date = cols[0]

        if "CREATED" in uprline:
            record.created = date, version
        elif "LAST SEQUENCE UPDATE" in uprline:
            record.sequence_update = date, version
        elif "LAST ANNOTATION UPDATE" in uprline:
            record.annotation_update = date, version
        else:
            raise SwissProtParserError("Unrecognised DT (DaTe) line", line=line)
    elif (
        "INTEGRATED INTO" in uprline
        or "SEQUENCE VERSION" in uprline
        or "ENTRY VERSION" in uprline
    ):
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        
        
        try:
            version = 0
            for s in cols[-1].split("."):
                if s.isdigit():
                    version = int(s)
        except ValueError:
            version = 0
        date = cols[0].rstrip(",")

        
        
        if "INTEGRATED" in uprline:
            record.created = date, version
        elif "SEQUENCE VERSION" in uprline:
            record.sequence_update = date, version
        elif "ENTRY VERSION" in uprline:
            record.annotation_update = date, version
        else:
            raise SwissProtParserError("Unrecognised DT (DaTe) line", line=line)
    else:
        raise SwissProtParserError("Failed to parse DT (DaTe) line", line=line)


def _read_ox(record, line):
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    line = line.split("{")[0]
    if record.taxonomy_id:
        ids = line[5:].rstrip().rstrip(";")
    else:
        descr, ids = line[5:].rstrip().rstrip(";").split("=")
        assert descr == "NCBI_TaxID", f"Unexpected taxonomy type {descr}"
    record.taxonomy_id.extend(ids.split(", "))


def _read_oh(record, line):
    
    assert line[5:].startswith("NCBI_TaxID="), f"Unexpected {line}"
    line = line[16:].rstrip()
    assert line[-1] == "." and line.count(";") == 1, line
    taxid, name = line[:-1].split(";")
    record.host_taxonomy_id.append(taxid.strip())
    record.host_organism.append(name.strip())


def _read_rn(reference, rn):
    
    
    
    
    words = rn.split(None, 1)
    number = words[0]
    assert number.startswith("[") and number.endswith("]"), f"Missing brackets {number}"
    reference.number = int(number[1:-1])
    if len(words) > 1:
        evidence = words[1]
        assert evidence.startswith("{") and evidence.endswith(
            "}"
        ), f"Missing braces {evidence}"
        reference.evidence = evidence[1:-1].split("|")


def _read_rc(reference, value):
    cols = value.split(";")
    if value[-1] == ";":
        unread = ""
    else:
        cols, unread = cols[:-1], cols[-1]
    for col in cols:
        if not col:  
            return
        
        i = col.find("=")
        if i >= 0:
            token, text = col[:i], col[i + 1 :]
            comment = token.lstrip(), text
            reference.comments.append(comment)
        else:
            comment = reference.comments[-1]
            comment = f"{comment} {col}"
            reference.comments[-1] = comment
    return unread


def _read_rx(reference, value):
    
    
    

    
    
    
    
    value = value.replace(" [NCBI, ExPASy, Israel, Japan]", "")

    
    
    
    
    
    
    
    warn = False
    if "=" in value:
        cols = value.split("; ")
        cols = [x.strip() for x in cols]
        cols = [x for x in cols if x]
        for col in cols:
            x = col.split("=")
            if len(x) != 2 or x == ("DOI", "DOI"):
                warn = True
                break
            assert len(x) == 2, f"I don't understand RX line {value}"
            reference.references.append((x[0], x[1].rstrip(";")))
    
    else:
        cols = value.split("; ")
        
        if len(cols) != 2:
            warn = True
        else:
            reference.references.append((cols[0].rstrip(";"), cols[1].rstrip(".")))
    if warn:
        import warnings

        from Bio import BiopythonParserWarning

        warnings.warn(f"Possibly corrupt RX line {value!r}", BiopythonParserWarning)


def _read_cc(record, line):
    key, value = line[5:8], line[9:].rstrip()
    if key == "-!-":  
        record.comments.append(value)
    elif key == "   ":  
        if not record.comments:
            
            record.comments.append(value)
        else:
            record.comments[-1] += " " + value


def _read_dr(record, value):
    cols = value.rstrip(".").split("; ")
    record.cross_references.append(tuple(cols))


def _read_pe(record, value):
    pe = value.split(":")
    record.protein_existence = int(pe[0])


def _read_kw(record, value):
    
    
    
    
    
    
    
    
    
    
    for val in value.rstrip(";.").split("; "):
        if val.endswith("}"):
            
            val = val.rsplit("{", 1)[0]
        record.keywords.append(val.strip())


def _read_ft(record, line):
    name = line[5:13].rstrip()
    if name:
        if line[13:21] == "        ":  
            location = line[21:80].rstrip()
            try:
                isoform_id, location = location.split(":")
            except ValueError:
                isoform_id = None
            try:
                from_res, to_res = location.split("..")
            except ValueError:
                from_res = location
                to_res = ""
            qualifiers = {}
        else:  
            from_res = line[14:20].lstrip()
            to_res = line[21:27].lstrip()
            isoform_id = None
            description = line[34:75].rstrip()
            qualifiers = {"description": description}
        from_res = Position.fromstring(from_res, -1)
        if to_res == "":
            to_res = from_res + 1
        else:
            to_res = Position.fromstring(to_res)
        location = SimpleLocation(from_res, to_res, ref=isoform_id)
        feature = FeatureTable(
            location=location, type=name, id=None, qualifiers=qualifiers
        )
        record.features.append(feature)
        return
    
    feature = record.features[-1]
    if line[5:34] == "                             ":  
        description = line[34:75].rstrip()
        if description.startswith("/FTId="):
            
            feature.id = description[6:].rstrip(".")
            return
        
        old_description = feature.qualifiers["description"]
        if old_description.endswith("-"):
            description = f"{old_description}{description}"
        else:
            description = f"{old_description} {description}"

        if feature.type in ("VARSPLIC", "VAR_SEQ"):  
            
            
            
            
            
            try:
                first_seq, second_seq = description.split(" -> ")
            except ValueError:
                pass
            else:
                extra_info = ""
                
                
                extra_info_pos = second_seq.find(" (")
                if extra_info_pos != -1:
                    extra_info = second_seq[extra_info_pos:]
                    second_seq = second_seq[:extra_info_pos]
                
                first_seq = first_seq.replace(" ", "")
                second_seq = second_seq.replace(" ", "")
                
                description = first_seq + " -> " + second_seq + extra_info
        feature.qualifiers["description"] = description
    else:  
        value = line[21:].rstrip()
        match = re.match(r"^/([a-z_]+)=", value)
        if match:
            qualifier_type = match.group(1)
            value = value[len(match.group(0)) :]
            if not value.startswith('"'):
                raise ValueError("Missing starting quote in feature")
            if qualifier_type == "id":
                if not value.endswith('"'):
                    raise ValueError("Missing closing quote for id")
                feature.id = value[1:-1]
            else:
                if value.endswith('"'):
                    value = value[1:-1]
                else:  
                    value = value[1:]
                if qualifier_type in feature.qualifiers:
                    raise ValueError(
                        f"Feature qualifier {qualifier_type!r} already exists for feature"
                    )
                feature.qualifiers[qualifier_type] = value
            return
        
        keys = list(feature.qualifiers.keys())
        key = keys[-1]
        description = value.rstrip('"')
        old_description = feature.qualifiers[key]
        if key == "evidence" or old_description.endswith("-"):
            description = f"{old_description}{description}"
        else:
            description = f"{old_description} {description}"
        if feature.type == "VAR_SEQ":  
            try:
                first_seq, second_seq = description.split(" -> ")
            except ValueError:
                pass
            else:
                extra_info = ""
                
                
                extra_info_pos = second_seq.find(" (")
                if extra_info_pos != -1:
                    extra_info = second_seq[extra_info_pos:]
                    second_seq = second_seq[:extra_info_pos]
                
                first_seq = first_seq.replace(" ", "")
                second_seq = second_seq.replace(" ", "")
                
                description = first_seq + " -> " + second_seq + extra_info
        feature.qualifiers[key] = description


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
