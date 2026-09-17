





"""Bio.Align support for "clustal" output from CLUSTAL W and other tools.

You are expected to use this module via the Bio.Align functions (or the
Bio.SeqIO functions if you are interested in the sequences only).
"""

import Bio
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Clustalw alignment writer."""

    fmt = "Clustal"

    def write_header(self, stream, alignments):
        """Use this to write the file header."""
        try:
            metadata = alignments.metadata
            program = metadata["Program"]
        except (AttributeError, KeyError):
            program = "Biopython"
            version = Bio.__version__
        else:
            version = metadata.get("Version", "")
        line = f"{program} {version} multiple sequence alignment\n"
        stream.write(line)
        stream.write("\n")
        stream.write("\n")

    def format_alignment(self, alignment):
        """Return a string with a single alignment in the Clustal format."""
        nseqs, length = alignment.shape
        if nseqs == 0:
            raise ValueError("Must have at least one sequence")
        if length == 0:
            raise ValueError("Non-empty sequences are required")

        try:
            column_annotations = alignment.column_annotations
        except AttributeError:
            consensus = None
        else:
            consensus = column_annotations.get("clustal_consensus")

        gapped_sequences = list(alignment)
        names = []
        for i, sequence in enumerate(alignment.sequences):
            try:
                name = sequence.id
            except AttributeError:
                name = "sequence_%d" % i  
            else:
                
                
                
                
                name = name[:30].replace(" ", "_")
            name = name.ljust(36)
            names.append(name)

        lines = []
        start = 0
        while start != length:
            
            
            stop = start + 50
            if stop > length:
                stop = length

            for name, gapped_sequence in zip(names, gapped_sequences):
                line = f"{name}{gapped_sequence[start:stop]}\n"
                lines.append(line)

            
            if consensus is not None:
                line = " " * 36 + consensus[start:stop] + "\n"
                lines.append(line)

            lines.append("\n")
            start = stop
        lines.append("\n")
        return "".join(lines)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Clustalw alignment iterator."""

    fmt = "Clustal"

    def _read_header(self, stream):
        line = stream.readline()
        if not line:
            raise ValueError("Empty file.")

        self.metadata = {}
        
        words = line.split()
        known_programs = [
            "CLUSTAL",
            "PROBCONS",
            "MUSCLE",
            "MSAPROBS",
            "Kalign",
            "Biopython",
        ]
        program = words[0]
        if program not in known_programs:
            raise ValueError(
                "%s is not known to generate CLUSTAL files: %s"
                % (program, ", ".join(known_programs))
            )
        self.metadata["Program"] = program

        
        for word in words:
            if word[0] == "(" and word[-1] == ")":
                word = word[1:-1]
            if word[0].isdigit():
                self.metadata["Version"] = word
                break

    def _read_next_alignment(self, stream):
        
        
        
        ids = []
        aligned_seqs = []
        consensus = ""
        index = None  

        
        for line in stream:
            if line.startswith(" "):
                
                assert len(ids) > 0
                assert index is not None
                length = len(aligned_seq)  
                consensus = line[index : index + length]
                break
            elif line.strip():
                
                fields = line.split()

                
                
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError("Could not parse line:\n%s" % line)

                seqid, aligned_seq = fields[:2]
                ids.append(seqid)
                aligned_seqs.append(aligned_seq)

                
                if index is None:
                    index = line.find(aligned_seq, len(seqid))

                if len(fields) == 3:
                    
                    try:
                        count = int(fields[2])
                    except ValueError:
                        raise ValueError(
                            "Could not parse line, bad sequence count:\n%s" % line
                        ) from None
                    if len(aligned_seq) - aligned_seq.count("-") != count:
                        raise ValueError(
                            "Could not parse line, incorrect sequence count:\n%s" % line
                        ) from None
            else:
                
                if index:
                    break
        else:
            return

        assert index is not None

        
        length = len(aligned_seqs[0])
        for aligned_seq in aligned_seqs:
            assert len(aligned_seq) == length
        if consensus:
            assert len(consensus) == length

        n = len(aligned_seqs)
        i = 0
        
        for line in stream:
            if line.startswith(" "):  
                assert index is not None
                length = len(aligned_seq)
                consensus += line[index : index + length]
            elif not line.strip():  
                continue
            else:
                seqid = ids[i]
                
                fields = line.split()

                
                
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError("Could not parse line:\n%s" % line)

                assert seqid == fields[0]
                aligned_seq = fields[1]
                aligned_seqs[i] += aligned_seq

                if len(fields) == 3:
                    
                    try:
                        count = int(fields[2])
                    except ValueError:
                        raise ValueError(
                            "Could not parse line, bad sequence number:\n%s" % line
                        ) from None
                    if len(aligned_seqs[i]) - aligned_seqs[i].count("-") != count:
                        raise ValueError(
                            "Could not parse line, incorrect sequence count:\n%s" % line
                        ) from None
                i += 1
                if i == n:
                    i = 0
        aligned_seqs = [s.encode() for s in aligned_seqs]
        seqs, coordinates = Alignment.parse_printed_alignment(aligned_seqs)
        records = [
            SeqRecord(Seq(seq), id=seqid, description="")
            for (seqid, seq) in zip(ids, seqs)
        ]
        alignment = Alignment(records, coordinates)
        if consensus:
            columns = alignment.length
            if len(consensus) != columns:
                raise ValueError(
                    "Alignment has %i columns, consensus length is %i, '%s'"
                    % (columns, len(consensus), consensus)
                )
            alignment.column_annotations = {}
            alignment.column_annotations["clustal_consensus"] = consensus
        return alignment
