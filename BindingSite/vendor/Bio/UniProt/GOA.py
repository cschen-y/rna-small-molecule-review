







"""Parsers for the GAF, GPA and GPI formats from UniProt-GOA.

Uniprot-GOA README + GAF format description:
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README

Gene Association File, GAF formats:
http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
http://geneontology.org/docs/go-annotation-file-gaf-format-2.0/

Gene Product Association Data  (GPA format) README:
http://geneontology.org/docs/gene-product-association-data-gpad-format/

Gene Product Information (GPI format) README:
http://geneontology.org/docs/gene-product-information-gpi-format/

Go Annotation files are located here:
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/
"""

import copy





GAF20FIELDS = [
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "Evidence",
    "With",
    "Aspect",
    "DB_Object_Name",
    "Synonym",
    "DB_Object_Type",
    "Taxon_ID",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID",
]


GAF10FIELDS = [
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "Evidence",
    "With",
    "Aspect",
    "DB_Object_Name",
    "Synonym",
    "DB_Object_Type",
    "Taxon_ID",
    "Date",
    "Assigned_By",
]


GPA10FIELDS = [
    "DB",
    "DB_Object_ID",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "Evidence code",
    "With",
    "Interacting_taxon_ID",
    "Date",
    "Assigned_by",
    "Annotation_Extension",
    "Spliceform_ID",
]


GPA11FIELDS = [
    "DB",
    "DB_Object_ID",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "ECO_Evidence_code",
    "With",
    "Interacting_taxon_ID",
    "Date",
    "Assigned_by",
    "Annotation Extension",
    "Annotation_Properties",
]


GPI10FIELDS = [
    "DB",
    "DB_subset",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "DB_Object_Name",
    "DB_Object_Synonym",
    "DB_Object_Type",
    "Taxon",
    "Annotation_Target_Set",
    "Annotation_Completed",
    "Parent_Object_ID",
]


GPI11FIELDS = [
    "DB_Object_ID",
    "DB_Object_Symbol",
    "DB_Object_Name",
    "DB_Object_Synonym",
    "DB_Object_Type",
    "Taxon",
    "Parent_Object_ID",
    "DB_Xref",
    "Gene_Product_Properties",
]


GPI12FIELDS = [
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "DB_Object_Name",
    "DB_Object_Synonym",
    "DB_Object_Type",
    "Taxon",
    "Parent_Object_ID",
    "DB_Xref",
    "Gene_Product_Properties",
]


def _gpi10iterator(handle):
    """Read GPI 1.0 format files (PRIVATE).

    This iterator is used to read a gp_information.goa_uniprot
    file which is in the GPI 1.0 format.
    """
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[5] = inrec[5].split("|")  
        inrec[8] = inrec[8].split("|")  
        yield dict(zip(GPI10FIELDS, inrec))


def _gpi11iterator(handle):
    """Read GPI 1.1 format files (PRIVATE).

    This iterator is used to read a gp_information.goa_uniprot
    file which is in the GPI 1.1 format.
    """
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[2] = inrec[2].split("|")  
        inrec[3] = inrec[3].split("|")  
        inrec[7] = inrec[7].split("|")  
        inrec[8] = inrec[8].split("|")  
        yield dict(zip(GPI11FIELDS, inrec))


def _gpi12iterator(handle):
    """Read GPI 1.2 format files (PRIVATE).

    This iterator is used to read a gp_information.goa_uniprot
    file which is in the GPI 1.2 format.
    """
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split("|")  
        inrec[4] = inrec[4].split("|")  
        inrec[8] = inrec[8].split("|")  
        inrec[9] = inrec[9].split("|")  
        yield dict(zip(GPI12FIELDS, inrec))


def gpi_iterator(handle):
    """Read GPI format files.

    This function should be called to read a
    gp_information.goa_uniprot file. At the moment, there is
    only one format, but this may change, so
    this function is a placeholder a future wrapper.
    """
    inline = handle.readline()
    if inline.strip() == "!gpi-version: 1.2":
        return _gpi12iterator(handle)
    elif inline.strip() == "!gpi-version: 1.1":
        
        return _gpi11iterator(handle)
    elif inline.strip() == "!gpi-version: 1.0":
        
        return _gpi10iterator(handle)
    elif inline.strip() == "!gpi-version: 2.1":
        
        
        raise NotImplementedError("Sorry, parsing GPI version 2 not implemented yet.")
    else:
        raise ValueError(f"Unknown GPI version {inline}\n")


def _gpa10iterator(handle):
    """Read GPA 1.0 format files (PRIVATE).

    This iterator is used to read a gp_association.*
    file which is in the GPA 1.0 format. Do not call directly. Rather,
    use the gpaiterator function.
    """
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[2] = inrec[2].split("|")  
        inrec[4] = inrec[4].split("|")  
        inrec[6] = inrec[6].split("|")  
        inrec[10] = inrec[10].split("|")  
        yield dict(zip(GPA10FIELDS, inrec))


def _gpa11iterator(handle):
    """Read GPA 1.1 format files (PRIVATE).

    This iterator is used to read a gp_association.goa_uniprot
    file which is in the GPA 1.1 format. Do not call directly. Rather
    use the gpa_iterator function
    """
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[2] = inrec[2].split("|")  
        inrec[4] = inrec[4].split("|")  
        inrec[6] = inrec[6].split("|")  
        inrec[10] = inrec[10].split("|")  
        yield dict(zip(GPA11FIELDS, inrec))


def gpa_iterator(handle):
    """Read GPA format files.

    This function should be called to read a
    gene_association.goa_uniprot file. Reads the first record and
    returns a gpa 1.1 or a gpa 1.0 iterator as needed
    """
    inline = handle.readline()
    if inline.strip() == "!gpa-version: 1.1":
        
        return _gpa11iterator(handle)
    elif inline.strip() == "!gpa-version: 1.0":
        
        return _gpa10iterator(handle)
    else:
        raise ValueError(f"Unknown GPA version {inline}\n")


def _gaf20iterator(handle):
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split("|")  
        inrec[5] = inrec[5].split("|")  
        inrec[7] = inrec[7].split("|")  
        inrec[10] = inrec[10].split("|")  
        inrec[12] = inrec[12].split("|")  
        yield dict(zip(GAF20FIELDS, inrec))


def _gaf10iterator(handle):
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split("|")  
        inrec[5] = inrec[5].split("|")  
        inrec[7] = inrec[7].split("|")  
        inrec[10] = inrec[10].split("|")  
        inrec[12] = inrec[12].split("|")  
        yield dict(zip(GAF10FIELDS, inrec))


def _gaf10byproteiniterator(handle):
    cur_id = None
    id_rec_list = []
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split("|")  
        inrec[5] = inrec[5].split("|")  
        inrec[7] = inrec[7].split("|")  
        inrec[10] = inrec[10].split("|")  
        inrec[12] = inrec[12].split("|")  
        cur_rec = dict(zip(GAF10FIELDS, inrec))
        if cur_rec["DB_Object_ID"] != cur_id and cur_id:
            ret_list = copy.copy(id_rec_list)
            id_rec_list = [cur_rec]
            cur_id = cur_rec["DB_Object_ID"]
            yield ret_list
        else:
            cur_id = cur_rec["DB_Object_ID"]
            id_rec_list.append(cur_rec)


def _gaf20byproteiniterator(handle):
    cur_id = None
    id_rec_list = []
    for inline in handle:
        if inline[0] == "!":
            continue
        inrec = inline.rstrip("\n").split("\t")
        if len(inrec) == 1:
            continue
        inrec[3] = inrec[3].split("|")  
        inrec[5] = inrec[5].split("|")  
        inrec[7] = inrec[7].split("|")  
        inrec[10] = inrec[10].split("|")  
        inrec[12] = inrec[12].split("|")  
        cur_rec = dict(zip(GAF20FIELDS, inrec))
        if cur_rec["DB_Object_ID"] != cur_id and cur_id:
            ret_list = copy.copy(id_rec_list)
            id_rec_list = [cur_rec]
            cur_id = cur_rec["DB_Object_ID"]
            yield ret_list
        else:
            cur_id = cur_rec["DB_Object_ID"]
            id_rec_list.append(cur_rec)


def gafbyproteiniterator(handle):
    """Iterate over records in a gene association file.

    Returns a list of all consecutive records with the same DB_Object_ID
    This function should be called to read a
    gene_association.goa_uniprot file. Reads the first record and
    returns a gaf 2.0 or a gaf 1.0 iterator as needed
    2016-04-09: added GAF 2.1 iterator & fixed bug in iterator assignment
    In the meantime GAF 2.1 uses the GAF 2.0 iterator
    """
    inline = handle.readline()
    if inline.strip() == "!gaf-version: 2.0":
        
        return _gaf20byproteiniterator(handle)
    elif inline.strip() == "!gaf-version: 1.0":
        
        return _gaf10byproteiniterator(handle)
    elif inline.strip() == "!gaf-version: 2.1":
        
        
        return _gaf20byproteiniterator(handle)
    elif inline.strip() == "!gaf-version: 2.2":
        
        
        
        
        return _gaf20byproteiniterator(handle)
    else:
        raise ValueError(f"Unknown GAF version {inline}\n")


def gafiterator(handle):
    """Iterate over a GAF 1.0 or 2.x file.

    This function should be called to read a
    gene_association.goa_uniprot file. Reads the first record and
    returns a gaf 2.x or a gaf 1.0 iterator as needed

    Example: open, read, interat and filter results.

    Original data file has been trimmed to ~600 rows.

    Original source ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/YEAST/goa_yeast.gaf.gz

    >>> from Bio.UniProt.GOA import gafiterator, record_has
    >>> Evidence = {'Evidence': set(['ND'])}
    >>> Synonym = {'Synonym': set(['YA19A_YEAST', 'YAL019W-A'])}
    >>> Taxon_ID = {'Taxon_ID': set(['taxon:559292'])}
    >>> with open('UniProt/goa_yeast.gaf', 'r') as handle:
    ...     for rec in gafiterator(handle):
    ...         if record_has(rec, Taxon_ID) and record_has(rec, Evidence) and record_has(rec, Synonym):
    ...             for key in ('DB_Object_Name', 'Evidence', 'Synonym', 'Taxon_ID'):
    ...                 print(rec[key])
    ...
    Putative uncharacterized protein YAL019W-A
    ND
    ['YA19A_YEAST', 'YAL019W-A']
    ['taxon:559292']
    Putative uncharacterized protein YAL019W-A
    ND
    ['YA19A_YEAST', 'YAL019W-A']
    ['taxon:559292']
    Putative uncharacterized protein YAL019W-A
    ND
    ['YA19A_YEAST', 'YAL019W-A']
    ['taxon:559292']

    """
    inline = handle.readline()
    if inline.strip() == "!gaf-version: 2.0":
        
        return _gaf20iterator(handle)
    elif inline.strip() == "!gaf-version: 2.1":
        
        
        return _gaf20iterator(handle)
    elif inline.strip() == "!gaf-version: 2.2":
        
        
        
        
        return _gaf20iterator(handle)
    elif inline.strip() == "!gaf-version: 1.0":
        
        return _gaf10iterator(handle)
    else:
        raise ValueError(f"Unknown GAF version {inline}\n")


def writerec(outrec, handle, fields=GAF20FIELDS):
    """Write a single UniProt-GOA record to an output stream.

    Caller should know the  format version. Default: gaf-2.0
    If header has a value, then it is assumed this is the first record,
    a header is written.
    """
    outstr = ""
    for field in fields[:-1]:
        if isinstance(outrec[field], list):
            for subfield in outrec[field]:
                outstr += subfield + "|"
            outstr = outstr[:-1] + "\t"
        else:
            outstr += outrec[field] + "\t"
    outstr += outrec[fields[-1]] + "\n"
    handle.write(outstr)


def writebyproteinrec(outprotrec, handle, fields=GAF20FIELDS):
    """Write a list of GAF records to an output stream.

    Caller should know the  format version. Default: gaf-2.0
    If header has a value, then it is assumed this is the first record,
    a header is written. Typically the list is the one read by fafbyproteinrec, which
    contains all consecutive lines with the same DB_Object_ID
    """
    for outrec in outprotrec:
        writerec(outrec, handle, fields=fields)


def record_has(inrec, fieldvals):
    """Accept a record, and a dictionary of field values.

    The format is {'field_name': set([val1, val2])}.
    If any field in the record has  a matching value, the function returns
    True. Otherwise, returns False.
    """
    retval = False
    for field in fieldvals:
        if isinstance(inrec[field], str):
            set1 = {inrec[field]}
        else:
            set1 = set(inrec[field])
        if set1 & fieldvals[field]:
            retval = True
            break
    return retval


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
