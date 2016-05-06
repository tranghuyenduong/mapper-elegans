import re

from formats.Alignment import AlignmentRecord

class BowtieRecord():

    def __init__(self, read, strand, ref, start_coord, mapping_qual, cigar, mate_ref, mate_offset, infer_insert_size, seq,
        read_qual, xa, md=None, nm=None):
        self.read = read
        self.read_count = int(re.search('.*-(\d+)$', self.read).group(1))

        self.original_strand = strand
        strand = int(strand)
        if(strand & 0x10 != 0):
            self.strand = '-'
        else:
            self.strand = '+'

        self.original_ref = ref
        self.ref = ref
        self.original_start_coord = start_coord
        self.start_coord = int(start_coord)
        self.mapping_qual = int(mapping_qual)
        self.cigar = cigar
        self.mate_ref = mate_ref
        self.mate_offset = mate_offset
        self.infer_insert_size = infer_insert_size
        self.seq = seq
        self.read_qual = read_qual
        self.end_coord = self.start_coord + len(self.seq)
        self.original_seq = seq
        self.xa = xa
        self.md = md
        self.nm = nm

    def to_alignment(self):
        return AlignmentRecord(*[
            self.ref,
            self.start_coord,
            self.end_coord,
            self.read,
            self.read_count,
            self.strand,
            self.original_seq,
            self
        ])

    @property
    def summary(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            self.read,
            self.original_strand,
            self.original_ref,
            self.original_start_coord,
            self.mapping_qual,
            self.cigar,
            self.mate_ref,
            self.mate_offset,
            self.infer_insert_size,
            self.original_seq,
            self.read_qual
            )
